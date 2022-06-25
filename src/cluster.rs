use ahash::AHashMap;
use clap::ArgEnum;

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
}

pub fn reorder(lhs: usize, rhs: usize) -> (usize, usize) {
    if lhs > rhs {
        (rhs, lhs)
    } else {
        (lhs, rhs)
    }
}
