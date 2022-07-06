use std::{path::Path, fs::File, io::{self, BufRead}};
use ahash::AHashMap;
use clap::ArgEnum;
use crate::state::AlnState;

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

#[derive(Debug)]
pub struct ClusteringResult {
    pub clusters: Vec<Vec<(u32, u32)>>,
}

impl ClusteringResult {
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

    pub fn mwt_am_score(&self, state: &AlnState, graph: &Graph) -> f64 {
        let mut score = 0.0;
        let mut cid : AHashMap<(u32, u32), usize> = AHashMap::new();
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

    pub fn from_plaintext<P>(path : P, graph : &Graph) -> anyhow::Result<Self> where P: AsRef<Path> {
        let file = File::open(path)?;
        let mut clusters : Vec<Vec<(u32, u32)>> = vec![];
        for l in io::BufReader::new(file).lines() {
            let line = l?;
            let mut cluster : Vec<(u32, u32)> = vec![];
            for tok in line.split(' ') {
                let node_id : usize = tok.parse()?;
                cluster.push(graph.node_pos[node_id]);
            }
            clusters.push(cluster);
        }
        Ok(Self {
            clusters,
        })
    }
}

pub fn reorder(lhs: usize, rhs: usize) -> (usize, usize) {
    if lhs > rhs {
        (rhs, lhs)
    } else {
        (lhs, rhs)
    }
}
