mod combined;
mod external;
// mod gcm;
mod aln;
mod cluster;
mod decompose;
mod exact_solver;
mod merge;
mod naive_upgma;
mod rt;
mod state;
mod utils;

use clap::{Parser, Subcommand};
use itertools::Itertools;
use naive_upgma::triplets_to_sims;
use ordered_float::NotNan;
use rand::Rng;
use std::path::PathBuf;
use tracing::{debug, error, info, span, warn, Level};

#[derive(Parser, Debug, Hash, PartialEq)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
}

#[derive(Subcommand, Debug, PartialEq, Hash)]
enum SubCommand {
    Merge {
        #[clap(short, long, multiple_values = true)]
        input: Vec<PathBuf>,
        #[clap(short, long, multiple_values = true)]
        glues: Vec<PathBuf>,
        #[clap(short, long, multiple_values = true)]
        weights: Vec<NotNan<f64>>,
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
    let mut g = cluster::Graph {
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

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    use std::time::Instant;
    let now = Instant::now();
    tracing_subscriber::fmt::init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()?;
    let args = Args::parse();
    match args.cmd {
        SubCommand::Merge {
            input,
            glues,
            weights,
            output,
        } => {
            let w = if weights.is_empty() {
                None
            } else {
                Some(weights)
            };
            info!("Analysis: merging alignments");
            info!(
                "Merging configuration (# alignments to merge, # glues, weights): {}, {}, {:?}",
                input.len(),
                glues.len(),
                w
            );
            combined::oneshot_merge_alignments(&input, &glues, &w, &output)
                .expect("Failed to merge alignments");
        }
    }
    info!("Total elapsed time: {:?}", now.elapsed());
    Ok(())
}
