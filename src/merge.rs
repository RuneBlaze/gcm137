use ahash::AHashMap;
use ordered_float::NotNan;

use crate::{aln::AlnProcessor, external::request_alignment};
use futures::executor::block_on;
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use seq_io::{fasta::Reader, BaseRecord};
use std::{
    collections::BTreeSet,
    error::Error,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Arc,
};
use tokio::{sync::Semaphore, task};

use crate::{
    naive_upgma::{ClusteringResult, Graph},
    state::AlnState,
};

pub fn state_from_constraints(constraint_alns: &[PathBuf]) -> anyhow::Result<AlnState> {
    let mut p = StateFromConstraints::default();
    for aln in constraint_alns {
        // cid : constraint id
        let mut reader = Reader::from_path(aln)?;
        while let Some(result) = reader.next() {
            let rec = result?;
            p.on_record(&rec)?;
        }
        p.next_aln();
    }
    return Ok(p.take());
}

pub struct StateFromConstraints {
    state: AlnState,
    sequence_id: usize,
    columns: usize,
    cid: usize,
}

impl StateFromConstraints {
    pub fn next_aln(&mut self) {
        self.state.column_counts.push(self.columns);
        self.columns = 0;
        self.cid += 1;
    }
}

impl Default for StateFromConstraints {
    fn default() -> Self {
        Self {
            state: AlnState::new(),
            sequence_id: 0,
            columns: 0,
            cid: 0,
        }
    }
}

impl AlnProcessor for StateFromConstraints {
    type Output = AlnState;

    fn on_record(&mut self, rec: &seq_io::fasta::RefRecord) -> anyhow::Result<()> {
        let name = String::from_utf8(rec.head().iter().copied().collect_vec())?;
        let mut column = 0usize;
        self.state.names.push(name.clone());
        self.state.names2id.insert(name, self.sequence_id);
        // res.id2constraint
        let mut s_slice: Vec<(u32, u32)> = vec![];
        for l in rec.seq_lines() {
            for &c in l {
                if c != b'-' {
                    s_slice.push((self.cid as u32, column as u32));
                }
                column += 1;
            }
        }
        if self.columns <= 0 {
            self.columns = column;
        } else {
            assert_eq!(column, self.columns);
        }
        self.state.s.push(s_slice);
        self.sequence_id += 1;
        Ok(())
    }

    fn take(&mut self) -> Self::Output {
        std::mem::replace(&mut self.state, AlnState::new())
    }
}

type SparseGraph = AHashMap<(u32, u32), AHashMap<(u32, u32), f64>>;

pub fn build_subgraph(state: &AlnState, glue: &PathBuf) -> anyhow::Result<SparseGraph> {
    let s = &state.s;
    let mut res = AHashMap::default();
    let mut colors: Vec<AHashMap<(u32, u32), usize>> = vec![];
    let mut reader = Reader::from_path(glue)?;
    while let Some(result) = reader.next() {
        let rec = result?;
        let name = String::from_utf8(rec.head().iter().copied().collect_vec())?;
        let mut column = 0;
        let mut first_ele = false;
        if colors.is_empty() {
            first_ele = true;
        }
        let mut non_gap = 0;
        for l in rec.seq_lines() {
            for &c in l {
                if first_ele {
                    colors.push(AHashMap::default());
                }
                if c != b'-' {
                    let id = state.names2id[&name];
                    let c = s[id][non_gap];
                    let entry = colors[column].entry(c).or_default();
                    *entry += 1;
                    non_gap += 1;
                }
                column += 1;
            }
        }
    }
    for c in &colors {
        for (c1, c2) in c.keys().tuple_combinations() {
            let (oc1, oc2) = if *c1 > *c2 { (c2, c1) } else { (c1, c2) };
            let entry = res.entry(*oc1).or_insert_with(AHashMap::default);
            let entry2 = entry.entry(*oc2).or_default();
            *entry2 += (c[c1] * c[c2]) as f64;
        }
    }
    Ok(res)
}

pub fn build_graph(
    state: &AlnState,
    glues: &[PathBuf],
    weights: &Option<Vec<NotNan<f64>>>,
) -> anyhow::Result<Graph> {
    let subgraphs_: anyhow::Result<Vec<SparseGraph>> = glues
        .par_iter()
        .map(|glue| build_subgraph(state, glue))
        .collect();
    let subgraphs = subgraphs_?;
    // now we need to merge the subgraphs
    let mut merged: AHashMap<usize, AHashMap<usize, f64>> = AHashMap::default();
    let mut pos2id: AHashMap<(u32, u32), usize> = AHashMap::default();
    let mut labels: BTreeSet<usize> = BTreeSet::default();
    let mut node_pos: Vec<(u32, u32)> = vec![];
    let mut id = 0;
    for (c, &l) in state.column_counts.iter().enumerate() {
        for i in 0..l {
            pos2id.insert((c as u32, i as u32), id);
            node_pos.push((c as u32, i as u32));
            id += 1;
        }
    }
    for (i, subgraph) in subgraphs.iter().enumerate() {
        let subgraph_weight = weights.as_ref().map(|w| w[i as usize].into_inner()).unwrap_or(1.0);
        for (u, map) in subgraph {
            for (v, w) in map {
                let u = pos2id[u];
                let v = pos2id[v];
                labels.insert(u);
                labels.insert(v);
                let entry = merged.entry(u).or_insert_with(AHashMap::default);
                let entry2 = entry.entry(v).or_default();
                *entry2 += w * subgraph_weight;
            }
        }
    }
    Ok(Graph {
        size: id,
        labels: labels.into_iter().collect_vec(),
        sims: merged,
        node_pos,
    })
}

pub fn build_frames(state: &AlnState, res: &ClusteringResult) -> Vec<Vec<u32>> {
    let k = state.column_counts.len();
    let mut last_frontier = vec![-1i64; k];
    let mut skip_frames: Vec<Vec<u32>> = vec![Vec::new(); k];
    for i in 0..k {
        skip_frames[i].push(0);
    }
    let traces = &res.clusters;
    for tr in traces {
        let mut c = 0u32;
        for e in tr.iter().chain(&[(k as u32, 123123 as u32)]) {
            while c < e.0 && c < k as u32 {
                // c not present in tr
                let l = skip_frames[c as usize].len();
                skip_frames[c as usize][l - 1] += 1;
                c += 1;
            }
            if c >= k as u32 {
                break;
            }
            assert_eq!(c, e.0);
            // check for singletons
            if last_frontier[c as usize] < (e.1 as i64 - 1) {
                // then singletons exist, we need to take care of them
                let num_singletons = (e.1 as i64 - last_frontier[c as usize] - 1) as usize;
                for j in 0..k {
                    if j == c as usize {
                        for _ in 0..num_singletons {
                            skip_frames[c as usize].push(0);
                        }
                    } else {
                        let l = skip_frames[j as usize].len();
                        skip_frames[j][l - 1] += num_singletons as u32;
                    }
                }
            }
            assert!(
                last_frontier[c as usize] < e.1 as i64,
                "last_frontier[{}] = {} >= {}",
                c,
                last_frontier[c as usize],
                e.1
            );
            last_frontier[c as usize] = e.1 as i64;

            skip_frames[c as usize].push(0);
            c += 1;
        }
    }

    for c in 0..k {
        let expected_len = state.column_counts[c] - 1;
        if (last_frontier[c as usize] as usize) < expected_len {
            // then singletons exist, we need to take care of them
            let num_singletons = (expected_len - last_frontier[c as usize] as usize) as usize;
            for j in 0..k {
                if j == c as usize {
                    for _ in 0..num_singletons {
                        skip_frames[c as usize].push(0);
                    }
                } else {
                    let l = skip_frames[j as usize].len();
                    skip_frames[j][l - 1] += num_singletons as u32;
                }
            }
        }
    }
    return skip_frames;
}

pub fn merge_alignments_from_frames(
    constraints: &[PathBuf],
    frames: &[Vec<u32>],
    outfile: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let out = File::create(outfile).unwrap();
    let mut writer = BufWriter::new(out);
    for (constraint, frame) in constraints.iter().zip(frames) {
        let mut reader = Reader::from_path(constraint)?;
        while let Some(result) = reader.next() {
            let rec = result?;
            let mut buf: Vec<u8> = vec![];
            let mut char_count = 0;
            let mut written = 0;
            for l in rec.seq_lines() {
                for c in l {
                    for _ in 0..frame[char_count] as usize {
                        buf.push(b'-');
                    }
                    written += frame[char_count] as usize + 1;
                    buf.push(*c);
                    char_count += 1;
                }
            }
            if char_count < frame.len() {
                for _ in 0..(frame[frame.len() - 1] + 1) {
                    buf.push(b'-');
                }
            }
            // println!("buflen: {}, sum of frame: {}, {}/{}, {}", buf.len(), written, char_count, frame.len(), frame[frame.len() - 1]);
            writer.write_all(b">")?;
            writer.write_all(rec.head()).unwrap();
            writer.write_all(b"\n").unwrap();
            buf.chunks(60).for_each(|chunk| {
                writer.write_all(chunk).unwrap();
                writer.write_all(b"\n").unwrap();
            });
        }
    }
    Ok(())
}

pub async fn align_glues(glues: &[PathBuf], tokens: usize) -> Vec<PathBuf> {
    let mut join_handles = vec![];
    let semaphore = Arc::new(Semaphore::new(tokens));
    let mut outputs = vec![];
    for g in glues {
        let cg = g.clone();
        let permit = semaphore.clone().acquire_owned().await.unwrap();
        let mut outpath = g.clone();
        outpath.set_extension("aln");
        let opc = outpath.clone();
        join_handles.push(task::spawn_blocking(move || {
            request_alignment(&cg, &outpath).expect("failed to align");
            drop(permit);
        }));
        outputs.push(opc);
    }
    futures::future::join_all(join_handles).await;
    outputs
}
