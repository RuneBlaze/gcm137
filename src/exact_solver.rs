use std::collections::VecDeque;

use ahash::AHashMap;
use fixedbitset::FixedBitSet;
use ndarray::{Array, ShapeBuilder};
use tracing::{debug, error, info, span, warn, Level};

use crate::{
    cluster::{reorder, ClusteringResult, Graph},
    state::{self, AlnState},
};

/// a DP algorithm for solving the exact two row MWT-AM problem
pub fn solve_twocase_mwt(
    graph: &AHashMap<(u32, u32), AHashMap<(u32, u32), f64>>,
) -> ClusteringResult {
    let mut edges: Vec<((u32, u32), (u32, u32), f64)> = vec![];
    for (u, map) in graph {
        for (v, f) in map {
            if u.1 < v.1 {
                // left to right
                edges.push((*u, *v, *f));
            } else {
                edges.push((*v, *u, *f));
            }
        }
    }
    // sort edges by the left coordinate
    edges.sort_by(|a, b| a.1.cmp(&b.1));
    // println!("{:?}", edges);
    let mut mem: AHashMap<FixedBitSet, f64> = AHashMap::default();
    let initial_state = FixedBitSet::with_capacity(edges.len());
    let initial_boundary = (0u32, 0u32);
    let mut back = vec![0usize; edges.len()];
    let res = two_case_mwt(
        &mut mem,
        &edges,
        &mut back,
        &(initial_state, initial_boundary, 0),
    );
    // println!("{:?}", res);
    let mut pt = back[0];
    let mut trace: Vec<Vec<(u32, u32)>> = vec![];
    while pt < edges.len() && back[pt] > pt {
        let mut e = (edges[pt].0, edges[pt].1);
        e = if e.0 .0 < e.1 .0 { e } else { (e.1, e.0) };
        trace.push(vec![e.0, e.1]);
        pt = back[pt + 1];
    }
    // println!("{:?}", back);
    let cr = ClusteringResult { clusters: trace };
    cr.check_validity();
    cr
}

pub fn weight_from_graph(graph: &Graph, i: usize, j: usize) -> f64 {
    let h1 = graph.sims.get(&i);
    if let Some(h2) = h1 {
        if let Some(f) = h2.get(&j) {
            return *f;
        }
    }
    0.0
}

pub fn sw_algorithm(graph: &Graph, state: &AlnState) -> ClusteringResult {
    let n = state.column_counts[0];
    let m = state.column_counts[1];
    let mut s = Array::<f64, _>::zeros((n + 1, m + 1).f());
    let mut back = Array::<u8, _>::zeros((n + 1, m + 1).f());
    for i in 0..(n + 1) {
        for j in 0..(m + 1) {
            if i == 0 || j == 0 {
                s[[i, j]] = 0.0;
                continue;
            }
            let mut max = 0.0;
            let mut max_pt = 0u8;
            let values = [
                s[[i - 1, j - 1]] + weight_from_graph(graph, i - 1, n + j - 1),
                s[[i - 1, j]],
                s[[i, j - 1]],
            ];
            for (i, &v) in values.iter().enumerate() {
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
    // println!("SW {}", score);
    let mut matches: Vec<Vec<(u32, u32)>> = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        let pt = back[[i, j]];
        if pt == 0 {
            i -= 1;
            j -= 1;
            matches.push(vec![(0, i as u32), (1, j as u32)]);
        } else if pt == 1 {
            i -= 1;
        } else if pt == 2 {
            j -= 1;
        }
    }
    matches.reverse();
    // println!("{:?}", matches);
    ClusteringResult { clusters: matches }
}

fn can_take(boundary: (u32, u32), edge_x: (u32, u32)) -> bool {
    if edge_x.0 >= boundary.0 && edge_x.1 >= boundary.1 {
        return true;
    }
    return false;
}
// #[tracing::instrument]
fn two_case_mwt(
    mem: &mut AHashMap<FixedBitSet, f64>,
    edges: &Vec<((u32, u32), (u32, u32), f64)>,
    back: &mut [usize],
    state: &(FixedBitSet, (u32, u32), usize),
) -> f64 {
    let (taken, boundary, frontier) = state;
    let b = *boundary;
    // println!("{:?}", b);
    if let Some(res) = mem.get(taken) {
        return *res;
    }
    if frontier >= &edges.len() {
        return 0.0;
    }
    // println!("{:?} {:?}", taken, edges.len() -  frontier);
    let mut maximum: f64 = 0f64;
    let mut best_next = 0usize;
    for i in *frontier..(edges.len()) {
        let e = edges[i];
        let e_x: (u32, u32) = (e.0 .1, e.1 .1);
        if can_take(b, e_x) {
            let mut new_state = taken.clone();
            new_state.insert(i);
            new_state.set_range(0..i, false);
            let new_boundary = (e_x.0 + 1, e_x.1 + 1);
            let new_frontier = i + 1;
            let new_res =
                e.2 + two_case_mwt(mem, edges, back, &(new_state, new_boundary, new_frontier));
            if new_res > maximum {
                maximum = new_res;
                best_next = i;
            }
        }
    }
    back[*frontier] = best_next;
    mem.insert(taken.clone(), maximum);
    return maximum;
}
