use std::{sync::Arc};

use ahash::AHashMap;
use crossbeam::{sync::WaitGroup, thread};
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use ndarray::{Array, ShapeBuilder};
use parking_lot::{RwLock, Mutex};
use rand::seq::SliceRandom;
use rand::{prelude::SmallRng, Rng, SeedableRng};
use tracing::{info, span, Level, debug};

use crate::{
    cluster::{ClusteringResult, Graph},
    state::AlnState,
};

fn random_partition(rng: &mut SmallRng, k: usize, bitset: &mut FixedBitSet) {
    let lb = 1;
    let ub = k - lb;
    let partition_size = rng.gen_range(lb..ub);
    let ix = rand::seq::index::sample(rng, k, partition_size);
    for i in ix {
        bitset.set(i, true);
    }
}

pub fn iterative_refinement(state: &AlnState, graph: &Graph, mut res: ClusteringResult) -> ClusteringResult {
    
    let k = state.column_counts.len();
    // let mut rng = SmallRng::from_entropy();
    let mut partition = FixedBitSet::with_capacity(k);
    let og_score = res.mwt_am_score(state, graph) as u32;
    info!(starting_score = og_score, "iterative refinement started");
    let mut rest_its = 1000usize;
    let mut sol = Arc::new(Mutex::new(res.clone()));
    {
        let mut w = sol.lock();
        w.mwt_am = og_score;
    }
    let mut sol2 = Arc::new(Mutex::new(res));
    {
        let mut w = sol.lock();
        w.mwt_am = og_score;
    }
    let wg = WaitGroup::new();
    let each_group = 8usize;
    for g in 0..2 {
        thread::scope(|scope| {
            for i in 0..each_group {
                let span = span!(Level::INFO, "opt", tid = i, group = g);
                let wg = wg.clone();
                let s = if g == 0 { (&sol).clone() } else { (&sol2).clone() };
                
                scope.spawn(move |_| {
                    let _enter = span.enter();
                    let mut rng = SmallRng::from_entropy();
                    let mut part = FixedBitSet::with_capacity(k);
                    let mut frustration = -1 as i32;
                    for epoch in 0..100 {
                        let mut r = {
                            let read = s.lock();
                            read.clone()
                        };
                        if frustration == 3 || frustration == -1 {
                            let mut slice = (0..k).collect_vec();
                            let ub = match frustration {
                                3 => 1,
                                _ => 2,
                            };
                            for _ in 0..ub {
                                slice.shuffle(&mut rng);
                                for &i in &slice {
                                    part.set(i, true);
                                    iterative_refinement_step(state, graph, &mut r, &part);
                                    part.set(i, false);
                                }
                            }
                        } else {
                            for _ in 0..50 {
                                random_partition(&mut rng, k, &mut part);
                                iterative_refinement_step(state, graph, &mut r, &part);
                                part.clear();
                            }
                        }
                        let read = s.lock();
                        let score = r.mwt_am_score(state, graph) as u32;
                        let their = {
                            read.mwt_am
                        };
                        if score > their {
                            let mut w = s.lock();
                            w.mwt_am = score as u32;
                            info!(new_score = w.mwt_am, improve = format!("{:.3}%", (score as f64 / og_score as f64 - 1.0) * 100.0) , frustration, "score improved");
                            w.clusters = r.clusters;
                            frustration = 0;
                        } else {
                            frustration += 1;
                            debug!("did not find better solution, frustration: {}", frustration);
                        }
                    }
                    drop(wg);
                });
            }
        }).unwrap();
    }
    wg.wait();

    {
        let l1 = sol.lock();
        let l2 = sol2.lock();
        let mut w = if l1.mwt_am > l2.mwt_am { l1 } else { l2 };
        for c in &mut w.clusters {
            c.sort_unstable_by_key(|x| x.0);
        }
    }
    // let lock = *sol;
    Arc::try_unwrap(sol).unwrap().into_inner()
}

#[inline]
fn get_graph_sim(g: &AHashMap<usize, AHashMap<usize, f64>>, u: usize, v: usize) -> Option<f64> {
    g.get(&u).and_then(|m| m.get(&v)).copied()
}

fn iterative_refinement_step(
    state: &AlnState,
    graph: &Graph,
    res: &mut ClusteringResult,
    partition: &FixedBitSet,
) {
    let mut c1: Vec<Vec<(u32, u32)>> = vec![];
    let mut c2: Vec<Vec<(u32, u32)>> = vec![];
    let mut pos2cid: AHashMap<(u32, u32), usize> = AHashMap::default();

    for (_, tr) in res.clusters.iter().enumerate() {
        let mut c1_buf: Vec<(u32, u32)> = vec![];
        let mut c2_buf: Vec<(u32, u32)> = vec![];
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
    let n = c1.len();
    let m = c2.len();
    let mut sims = Array::<u32, _>::zeros((n, m).f());
    for (&u, map) in &graph.sims {
        for (&v, &value) in map {
            let p1 = graph.node_pos[u];
            let p2 = graph.node_pos[v];
            if partition[p1.0 as usize] != partition[p2.0 as usize] {
                let (i1, i2) = if partition[p1.0 as usize] {
                    (
                        match pos2cid.get(&p1) {
                            Some(i) => *i,
                            None => panic!("{} {}", p1.0, p1.1),
                        },
                        match pos2cid.get(&p2) {
                            Some(i) => *i,
                            None => panic!(
                                "{} {} {} -> {}, {}",
                                p2.0,
                                p2.1,
                                u,
                                v,
                                graph.labels.contains(&u) && graph.labels.contains(&v)
                            ),
                        },
                    )
                } else {
                    (
                        match pos2cid.get(&p2) {
                            Some(i) => *i,
                            None => panic!("{} {}", p2.0, p2.1),
                        },
                        match pos2cid.get(&p1) {
                            Some(i) => *i,
                            None => panic!("{} {}", p1.0, p1.1),
                        },
                    )
                };
                sims[[i1, i2]] += value as u32;
            }
        }
    }

    // now run Smith-Waterman
    let mut s = Array::<u32, _>::zeros((n + 1, m + 1).f());
    let mut back = Array::<u8, _>::zeros((n + 1, m + 1).f());
    for i in 0..(n + 1) {
        for j in 0..(m + 1) {
            if i == 0 || j == 0 {
                s[[i, j]] = 0;
                if i == 0 {
                    back[[i, j]] = 2;
                }
                if j == 0 {
                    back[[i, j]] = 1;
                }
                continue;
            }
            let mut max = 0;
            let mut max_pt = 0u8;
            let w = sims[[i - 1, j - 1]];
            let values = [s[[i - 1, j - 1]] + w, s[[i - 1, j]], s[[i, j - 1]]];
            for (i, &v) in values.iter().enumerate() {
                if i == 0 && w <= 0 {
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
    let score = s[[n, m]];
    // if score <= res.mwt_am {
    //     return None;
    // }
    let mut matches: Vec<Vec<(u32, u32)>> = Vec::new();
    let (mut i, mut j) = (n, m);
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
    // Some(ClusteringResult {
    //     clusters: matches,
    //     mwt_am: score,
    // })
}
