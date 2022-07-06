use ahash::AHashMap;
use fixedbitset::FixedBitSet;
use ndarray::{Array, ShapeBuilder};
use rand::Rng;

use crate::{
    cluster::{ClusteringResult, Graph},
    state::AlnState,
};

fn random_partition(k: usize, bitset: &mut FixedBitSet) {
    let lb = 1;
    let ub = k - lb;
    let mut rng = rand::thread_rng();
    let partition_size = rng.gen_range(lb..ub);
    let ix = rand::seq::index::sample(&mut rng, k, partition_size);
    for i in ix {
        bitset.set(i, true);
    }
}

pub fn iterative_refinement(state: &AlnState, graph: &Graph, res: &mut ClusteringResult) {
    let k = state.column_counts.len();
    let mut partition = FixedBitSet::with_capacity(k);
    let mut rest_its = 1000usize;
    for _ in 0..3 {
        for i in 0..k {
            partition.set(i, true);
            iterative_refinment_step(state, graph, res, &partition);
            partition.set(i, false);
        }
    }
    while rest_its > 0 {
        random_partition(k, &mut partition);
        iterative_refinment_step(state, graph, res, &partition);
        partition.clear();
        rest_its -= 1;
    }
    for c in &mut res.clusters {
        c.sort_unstable_by_key(|x| x.0);
    }
}

#[inline]
fn get_graph_sim(g : &AHashMap<usize, AHashMap<usize, f64>>, u : usize, v : usize) -> Option<f64> {
    g.get(&u).and_then(|m| m.get(&v)).copied()
}

fn iterative_refinment_step(
    state: &AlnState,
    graph: &Graph,
    res: &mut ClusteringResult,
    partition: &FixedBitSet,
) {
    let mut c1 : Vec<Vec<(u32, u32)>> = vec![];
    let mut c2 : Vec<Vec<(u32, u32)>> = vec![];
    let mut pos2cid : AHashMap<(u32, u32), usize> = AHashMap::default();
    let mut sims : AHashMap<usize, AHashMap<usize, f64>> = AHashMap::default();
    for (_, tr) in res.clusters.iter().enumerate() {
        let mut c1_buf : Vec<(u32, u32)> = vec![];
        let mut c2_buf : Vec<(u32, u32)> = vec![];
        for e in tr {
            if partition[e.0 as usize] {
                c1_buf.push(*e);
                pos2cid.insert((e.0, e.1), c1.len());
            } else {
                c2_buf.push(*e);
                pos2cid.insert((e.0, e.1), c2.len());
            }
        }
        if !c1_buf.is_empty() {
            c1.push(c1_buf);
        }
        if !c2_buf.is_empty() {
            c2.push(c2_buf);
        }
    }
    let mut old_score = 0.0;
    let mut upperbound = 0.0;
    for (&u, map) in &graph.sims {
        for (&v, &value) in map {
            let p1 = graph.node_pos[u];
            let p2 = graph.node_pos[v];
            if partition[p1.0 as usize] != partition[p2.0 as usize] {
                let (i1, i2) = if partition[p1.0 as usize] {
                    (match pos2cid.get(&p1) {
                        Some(i) => *i,
                        None => continue,
                    }, 
                    match pos2cid.get(&p2) {
                        Some(i) => *i,
                        None => continue,
                    })
                } else {
                    (match pos2cid.get(&p2) {
                        Some(i) => *i,
                        None => continue,
                    }, 
                    match pos2cid.get(&p1) {
                        Some(i) => *i,
                        None => continue,
                    })
                };
                let entry = sims.entry(i1).or_insert_with(AHashMap::default);
                let entry2 = entry.entry(i2).or_default();
                *entry2 += value;
                upperbound += value;
            }
        }
    }

    let n = c1.len();
    let m = c2.len();

    // now run Smith-Waterman
    let mut s = Array::<f64, _>::zeros((n + 1, m + 1).f());
    let mut back = Array::<u8, _>::zeros((n + 1, m + 1).f());
    for i in 0..(n + 1) {
        for j in 0..(m + 1) {
            if i == 0 || j == 0 {
                s[[i, j]] = 0.0;
                if i == 0 {
                    back[[i, j]] = 2;
                }
                if j == 0 {
                    back[[i, j]] = 1;
                }
                continue;
            }
            let mut max = 0.0;
            let mut max_pt = 0u8;
            let w = get_graph_sim(&sims, i - 1, j - 1).unwrap_or_default();
            let values = [s[[i - 1, j - 1]] + w, s[[i - 1, j]], s[[i, j - 1]]];
            for (i, &v) in values.iter().enumerate() {
                if i == 0 && w <= 0.0 {
                    max_pt = 1;
                    continue;
                }
                if v > max {
                    max = v;
                    max_pt = i as u8;
                }
            }
            s[[i, j]] = max;
            back[[i, j]] = max_pt;
        }
    }
    let _score = s[[n, m]];
    let mut matches: Vec<Vec<(u32, u32)>> = Vec::new();
    let (mut i, mut j) = (n, m);
    // let mut total_matches = 0;
    while !(i == 0 && j == 0) {
        let pt = back[[i, j]];
        let mut buf: Vec<(u32, u32)> = vec![];
        if pt == 0 {
            i -= 1;
            j -= 1;
            buf.extend(c1[i].iter().cloned());
            buf.extend(c2[j].iter().cloned());
        } else if pt == 1 {
            i -= 1;
            buf.extend(c1[i].iter().cloned());
        } else if pt == 2 {
            j -= 1;
            buf.extend(c2[j].iter().cloned());
        }
        matches.push(buf);
    }
    matches.reverse();
    res.clusters = matches;
}