use std::collections::{BTreeSet, BinaryHeap, VecDeque};

use crate::cluster::{ClusteringResult, Graph};
use crate::{cluster::reorder, state::AlnState};
use ahash::{AHashMap, AHashSet};
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use ordered_float::NotNan;
use petgraph::unionfind::UnionFind;
use roaring::RoaringBitmap;

pub fn triplets_to_sims(
    triplets: Vec<(usize, usize, f64)>,
) -> AHashMap<usize, AHashMap<usize, f64>> {
    let mut sims = AHashMap::new();
    for (i, j, sim) in triplets {
        sims.entry(i).or_insert(AHashMap::new());
        sims.entry(j).or_insert(AHashMap::new());
        sims.get_mut(&i).unwrap().insert(j, sim);
    }
    sims
}

pub fn union(lhs: &mut Vec<usize>, rhs: &[usize]) {
    let mut set = BTreeSet::from_iter(lhs.drain(0..));
    set.extend(rhs.iter());
    lhs.extend(set.into_iter());
}

fn ds_dfs(outedges: &Vec<Vec<usize>>, visited: &mut FixedBitSet, u: usize, v: usize) -> bool {
    let mut stack: Vec<usize> = vec![];
    for e in &outedges[u] {
        stack.push(*e);
    }
    while let Some(n) = stack.pop() {
        for &e in &outedges[n] {
            if !visited.contains(e) {
                if e == v {
                    return true;
                }
                stack.push(e);
                visited.insert(e);
            }
        }
    }
    false
}

fn order_clusters(absorbed: &FixedBitSet, outedges: &Vec<Vec<usize>>) -> Vec<usize> {
    let n = outedges.len();
    let mut indegs = vec![0; n];
    for (ix, i) in outedges.iter().enumerate() {
        if !absorbed.contains(ix) {
            for &j in i {
                indegs[j] += 1;
            }
        }
    }

    let mut stack: Vec<usize> = vec![];
    for (i, e) in indegs.iter().enumerate() {
        if !absorbed.contains(i) && *e == 0 {
            stack.push(i);
        }
    }

    let mut order: Vec<usize> = vec![];
    while let Some(t) = stack.pop() {
        order.push(t);
        for &j in &outedges[t] {
            indegs[j] -= 1;
            if indegs[j] == 0 {
                stack.push(j);
            }
        }
    }
    order
}

pub fn disjoint_set_po_exists(
    visited: &mut FixedBitSet,
    outedges: &Vec<Vec<usize>>,
    u: usize,
    v: usize,
    mid_visited: &FixedBitSet,
) -> bool {
    if outedges[u].contains(&v) && outedges[v].contains(&u) {
        return false;
    }
    visited.clear();

    match (mid_visited[u], mid_visited[v]) {
        (true, true) => {
            if ds_dfs(outedges, visited, u, v) || ds_dfs(outedges, visited, v, u) {
                return false;
            }
        }
        (false, false) => {
            if ds_dfs(outedges, visited, u, v) || ds_dfs(outedges, visited, v, u) {
                return false;
            }
        }
        (true, false) => {
            if ds_dfs(outedges, visited, v, u) {
                return false;
            }
        }
        (false, true) => {
            if ds_dfs(outedges, visited, u, v) {
                return false;
            }
        }
    }
    true
}

pub fn naive_upgma(graph: &Graph, state: &AlnState) -> ClusteringResult {
    let n = graph.labels.len();
    let m = graph.size;
    let mut clusters: UnionFind<usize> = UnionFind::new(n);
    let mut pq: BinaryHeap<(NotNan<f64>, usize, usize)> = BinaryHeap::new();
    let mut rows: Vec<RoaringBitmap> = vec![RoaringBitmap::new(); n];
    let mut cluster_sizes: Vec<usize> = vec![1; n];
    let mut node2init_cluster: Vec<usize> = vec![0; m];
    let mut weightmap: Vec<AHashMap<usize, f64>> = vec![AHashMap::new(); n];
    let mut visited = FixedBitSet::with_capacity(n);
    let mut mid_visited = FixedBitSet::with_capacity(n);
    for (i, l) in graph.labels.iter().enumerate() {
        node2init_cluster[*l] = i;
        rows[i].insert(graph.node_pos[*l].0);
    }
    for (&u, map) in &graph.sims {
        for (&v, &value) in map {
            if u == v {
                continue;
            }
            let mut lhs = node2init_cluster[u];
            let mut rhs = node2init_cluster[v];
            (lhs, rhs) = reorder(lhs, rhs);
            weightmap[lhs].entry(rhs).or_insert(value);
            weightmap[rhs].entry(lhs).or_insert(value);
            pq.push((NotNan::new(value).unwrap(), lhs, rhs));
        }
    }
    let mut order_outedges: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut order_inedges: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut bound = 0;
    let mut lengths = VecDeque::from_iter(state.column_counts.iter().copied());
    for i in 0..(n - 1) {
        let first_num = graph.labels[i];
        let second_num = graph.labels[i + 1];

        while first_num >= bound {
            bound += lengths.pop_front().unwrap();
        }

        if first_num < bound && second_num < bound {
            order_outedges[node2init_cluster[first_num]].push(node2init_cluster[second_num]);
            order_inedges[node2init_cluster[second_num]].push(node2init_cluster[first_num]);
        }
    }

    let mut absorbed = FixedBitSet::with_capacity(n);
    let mut invalidated: AHashSet<(usize, usize)> = AHashSet::new();

    while let Some((v, l_, r_)) = pq.pop() {
        let (l, r) = reorder(l_, r_);
        if absorbed[l] || absorbed[r] {
            continue;
        }

        if invalidated.contains(&(l, r)) || v != weightmap[l][&r] {
            continue;
        }

        if !rows[l].is_disjoint(&rows[r])
            || !disjoint_set_po_exists(&mut visited, &order_outedges, l, r, &mid_visited)
        {
            weightmap[l].remove(&r);
            weightmap[r].remove(&l);
            invalidated.insert((l, r));
            continue;
        }

        // println!("printing outedges...");
        // for (i, vec) in order_outedges.iter().enumerate() {
        //     for v in BTreeSet::from_iter(vec.iter()) {
        //         if !absorbed[i] && !absorbed[*v] {
        //             println!("oe {:?}[{}] -> {:?}[{}]", graph.node_pos[i],i, graph.node_pos[*v],*v);
        //         }
        //     }
        // }

        // let x = order_outedges.iter().map(|it| it.len()).sum::<usize>();
        // let y = order_inedges.iter().map(|it| it.len()).sum::<usize>();
        // assert_eq!(x, y);

        // println!("Joining {:?}", (graph.node_pos[l], graph.node_pos[r]));
        mid_visited = visited.clone();
        // visited.clear();
        clusters.union(l, r);
        let n = clusters.find(l);
        let m = if l == n { r } else { l };
        // println!("l={:?} r={:?} n={:?} m={:?}", l, r, n, m);
        absorbed.insert(m);
        // println!("BEFORE: oe: n {:?} oe: m{:?}", order_outedges[n], order_outedges[m]);
        // order_outedges.split_at()
        order_outedges.push(vec![]);
        let t1 = order_outedges.swap_remove(m);
        union(order_outedges[n].as_mut(), &t1);
        // println!("AFTER: oe: n {:?} oe: m{:?}", order_outedges[n], order_outedges[m]);
        // union(order_inedges[n].as_mut(), order_outedges[n].as_mut());
        // println!()

        // println!("BEFORE: ie: n {:?} ie: m{:?}", order_inedges[n], order_inedges[m]);
        order_inedges.push(vec![]);
        let t2 = order_inedges.swap_remove(m);
        union(order_inedges[n].as_mut(), &t2);
        // println!("AFTER: ie: n {:?} ie: m{:?}", order_inedges[n], order_inedges[m]);

        for &innode in &t2 {
            for e in order_outedges[innode].iter_mut() {
                if *e == m {
                    *e = n;
                }
            }
        }
        for &outnode in &t1 {
            for e in order_inedges[outnode].iter_mut() {
                if *e == m {
                    *e = n;
                }
            }
        }
        let mut seen: BTreeSet<usize> = BTreeSet::new();
        seen.extend(weightmap[l].keys().into_iter().chain(weightmap[r].keys()));
        for c in seen {
            if weightmap[l].contains_key(&c) && weightmap[r].contains_key(&c) {
                let v1 = weightmap[l][&c];
                let v2 = weightmap[r][&c];
                let pt = weightmap[n].get_mut(&c).unwrap();
                *pt = (v1 * cluster_sizes[l] as f64 + v2 * cluster_sizes[r] as f64)
                    / (cluster_sizes[l] + cluster_sizes[r]) as f64;
            } else if weightmap[l].contains_key(&c) {
                let v1 = weightmap[l][&c];
                weightmap[n].insert(c, v1);
            } else {
                let v2 = weightmap[r][&c];
                weightmap[n].insert(c, v2);
            }
            let v = weightmap[n][&c];
            weightmap[c].insert(n, v);
            let pair = reorder(c, n);
            pq.push((NotNan::new(v).unwrap(), pair.0, pair.1));
        }

        cluster_sizes[n] += cluster_sizes[m];
        rows.push(RoaringBitmap::new());
        let rm = rows.swap_remove(m);
        rows[n] |= rm;
    } // end UPGMA

    let cluster_ord = order_clusters(&absorbed, &order_outedges);
    let mut final_clusters: AHashMap<usize, Vec<usize>> = AHashMap::new();
    for i in &graph.labels {
        let cid = clusters.find(node2init_cluster[*i]);
        final_clusters.entry(cid).or_insert(Vec::new()).push(*i);
    }
    let mut ordered_clusters: Vec<Vec<(u32, u32)>> = Vec::new();
    for c in &cluster_ord {
        ordered_clusters.push(
            final_clusters[c]
                .iter()
                .map(|it| graph.node_pos[*it])
                .collect_vec(),
        );
    }

    // let mut score = 0.0;
    // for c in &cluster_ord {
    //     let x = &final_clusters[c];
    //     if x.len() <= 1 {
    //         continue;
    //     }
    //     let n1 = x[0];
    //     let n2 = x[1];
    //     score += graph.sims[&n1][&n2];
    // }
    // println!("score: {}", score);

    // for v in final_clusters.values() {

    // }
    ClusteringResult {
        clusters: ordered_clusters,
    }
}
