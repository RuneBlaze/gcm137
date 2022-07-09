use std::{collections::BTreeSet, path::PathBuf};

use fixedbitset::FixedBitSet;
use itertools::Itertools;
use seq_io::{fasta::Reader, BaseRecord};
use serde::Serialize;
use tracing::debug;

use crate::{
    cluster::ClusteringResult,
    merge::{load_graph, state_from_constraints},
    state::AlnState,
};

#[derive(Debug, Serialize)]
pub struct ScorerOutput {
    num_clusters: usize,
    mwt_am: f64,
}

impl ScorerOutput {
    pub fn new(num_clusters: usize, mwt_am: f64) -> Self {
        Self {
            num_clusters,
            mwt_am,
        }
    }
}

pub fn trace_from_alignment(
    state: &AlnState,
    alignment_path: &PathBuf,
) -> anyhow::Result<ClusteringResult> {
    let k = state.ncols();
    let mut trace: Vec<BTreeSet<(u32, u32)>> = vec![];
    let mut reader = Reader::from_path(alignment_path)?;
    let mut i = 0usize;
    let s = &state.s;

    while let Some(result) = reader.next() {
        let mut non_gap = 0;
        let rec = result?;
        let mut column = 0;
        let name = String::from_utf8(rec.head().iter().copied().collect_vec())?;
        let id = state.names2id[&name];
        for l in rec.seq_lines() {
            for &c in l {
                if i == 0 {
                    trace.push(BTreeSet::new());
                }
                if c != b'-' {
                    // println!("{} {} {}", i, column, std::str::from_utf8(&[c as u8])?);
                    let c = s.get(id).expect("sequence name not found!")[non_gap];
                    trace[column].insert(c);
                    non_gap += 1;
                }
                column += 1;
            }
        }
        i += 1;
    }
    let converted = trace
        .into_iter()
        .flat_map(|it| {
            if it.is_empty() {
                None
            } else {
                Some(it.into_iter().collect_vec())
            }
        })
        .collect_vec();
    let res = ClusteringResult::new(converted);
    Ok(res)
}

pub fn oneshot_score_alignment(
    constraints: &[PathBuf],
    graph_path: &PathBuf,
    alignment_path: &PathBuf,
) -> anyhow::Result<ScorerOutput> {
    let state = state_from_constraints(constraints)?;
    debug!("Constructed state from constraints");
    let graph = load_graph(&state, graph_path)?;
    debug!("Loaded graph");
    let trace = trace_from_alignment(&state, alignment_path)?;
    let score = trace.mwt_am_score(&state, &graph);
    Ok(ScorerOutput::new(trace.clusters.len(), score))
}
