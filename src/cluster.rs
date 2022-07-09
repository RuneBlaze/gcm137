use crate::state::AlnState;
use ahash::AHashMap;
use clap::ArgEnum;
use itertools::Itertools;
use std::{
    collections::BTreeSet,
    fs::File,
    io::{self, BufRead},
    mem::swap,
    path::Path,
};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug, Hash)]
pub enum GCMStep {
    Auto,
    Upgma,
    Pairwise,
}

pub struct Graph {
    pub size: usize,
    pub labels: Vec<usize>,
    pub sims: AHashMap<usize, AHashMap<usize, f64>>,
    pub node_pos: Vec<(u32, u32)>,
}

#[derive(Debug, Clone)]
pub struct ClusteringResult {
    pub clusters: Vec<Vec<(u32, u32)>>,
    pub mwt_am: u32,
}

impl ClusteringResult {
    pub fn new(sol: Vec<Vec<(u32, u32)>>) -> Self {
        ClusteringResult {
            clusters: sol,
            mwt_am: 0,
        }
    }
    pub fn check_validity(&self) {
        for i in 0..self.clusters.len() {
            for j in i..self.clusters.len() {
                if i == j {
                    continue;
                }
                let cond = Self::all_smaller(&self.clusters[i], &self.clusters[j])
                    || Self::all_smaller(&self.clusters[j], &self.clusters[i]);
                if !cond {
                    println!("{:?} {:?}", &self.clusters[i], &self.clusters[j]);
                    panic!("Clusters are not ordered");
                }
            }
        }
    }

    fn all_smaller(lhs: &[(u32, u32)], rhs: &[(u32, u32)]) -> bool {
        // checks if the x-coordinate of lhs is all smaller than the x-coordinate of rhs
        for (i, x) in lhs {
            if rhs.iter().any(|(j, y)| *i == *j && *x > *y) {
                return false;
            }
        }
        true
    }

    pub fn mwt_am_score(&self, _state: &AlnState, graph: &Graph) -> f64 {
        let mut score = 0.0;
        let mut cid: AHashMap<(u32, u32), usize> = AHashMap::new();
        for (i, tr) in self.clusters.iter().enumerate() {
            for e in tr {
                cid.insert(*e, i);
            }
        }
        for (&u, map) in &graph.sims {
            for (&v, &value) in map {
                let p1 = graph.node_pos[u];
                let p2 = graph.node_pos[v];
                if p1 == p2 {
                    continue;
                }
                let c1 = cid.get(&p1);
                let c2 = cid.get(&p2);
                match (c1, c2) {
                    (Some(i1), Some(i2)) if i1 == i2 => {
                        score += value;
                    }
                    (_, _) => continue,
                }
            }
        }
        score
    }

    /// re-insert singleton clusters that have weights into the trace
    pub fn hydrate(&mut self, state: &AlnState, graph: &Graph) {
        let mut lanes: Vec<BTreeSet<(u32, u32)>> = vec![BTreeSet::new(); state.ncols()];
        let labels = &graph.labels;
        for l in labels {
            let pos = graph.node_pos[*l];
            lanes[pos.0 as usize].insert(pos);
        }
        let mut lanes_order = lanes.iter().map(|it| it.iter().peekable()).collect_vec();
        let mut trace = vec![];
        swap(&mut trace, &mut self.clusters);
        for tr in trace {
            for e in &tr {
                // let mut it = &;
                while let Some(x) = lanes_order[e.0 as usize].next_if(|x| x.1 <= e.1) {
                    if x.1 < e.1 {
                        self.clusters.push(vec![*x]);
                    } else {
                        // info!("{:?}", x);
                    }
                }
            }
            self.clusters.push(tr);
        }
        for k in 0..state.ncols() {
            for x in lanes_order[k].by_ref() {
                self.clusters.push(vec![*x]);
            }
        }
    }

    pub fn from_plaintext<P>(path: P, graph: &Graph) -> anyhow::Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        let mut clusters: Vec<Vec<(u32, u32)>> = vec![];
        for l in io::BufReader::new(file).lines() {
            let line = l?;
            let mut cluster: Vec<(u32, u32)> = vec![];
            for tok in line.split(' ') {
                let node_id: usize = tok.parse()?;
                cluster.push(graph.node_pos[node_id]);
            }
            clusters.push(cluster);
        }
        Ok(Self::new(clusters))
    }
}

pub fn reorder(lhs: usize, rhs: usize) -> (usize, usize) {
    if lhs > rhs {
        (rhs, lhs)
    } else {
        (lhs, rhs)
    }
}
