mod combined;
mod external;
mod gcm;
mod merge;
mod naive_upgma;
mod state;
use clap::{Parser, Subcommand};
use itertools::Itertools;
use naive_upgma::triplets_to_sims;
use rand::Rng;
use std::path::PathBuf;

#[derive(Parser, Debug, Hash, PartialEq)]
#[clap(author, version, about)]
struct Args {
    // #[clap(short, long)]
    // input: PathBuf,
    // #[clap(short, long)]
    // tree: PathBuf,
    // #[clap(short, long)]
    // output: Option<PathBuf>,
    #[clap(subcommand)]
    cmd: SubCommand,
}

#[derive(Subcommand, Debug, Hash, PartialEq)]
enum SubCommand {
    Merge {
        #[clap(short, long, multiple_values = true)]
        input: Vec<PathBuf>,
        #[clap(short, long, multiple_values = true)]
        glues: Vec<PathBuf>,
        #[clap(short, long)]
        output: PathBuf,
    },
}

fn random_graph() {
    let lengths = [600; 25];
    let state = state::AlnState {
        names: vec![],
        names2id: ahash::AHashMap::default(),
        s: vec![],
        column_counts: Vec::from_iter(lengths.iter().copied()),
    };
    let size: usize = lengths.iter().sum();
    let mut g = naive_upgma::Graph {
        size,
        labels: (0..size).collect_vec(),
        sims: ahash::AHashMap::default(),
        node_pos: vec![],
    };
    for (i, l) in lengths.iter().enumerate() {
        for c in 0..*l {
            g.node_pos.push((i as u32, c as u32));
        }
    }
    let mut rng = rand::thread_rng();
    let mut triplets = vec![];
    for _ in 1..500000 {
        let x = rng.gen_range(0..size);
        let y = rng.gen_range(0..size);
        let f = rng.gen_range(1..5usize) as f64;
        triplets.push((x, y, f));
    }
    g.sims = triplets_to_sims(triplets);
    let res = naive_upgma::naive_upgma(&g, &state);
    res.check_validity();
    println!("{:?}", res);
}

fn main() {
    let args = Args::parse();
    match args.cmd {
        SubCommand::Merge {
            input,
            glues,
            output,
        } => {
            combined::oneshot_merge_alignments(&input, &glues, &output)
                .expect("Failed to merge alignments");
        }
    }
}