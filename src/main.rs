mod aln;
mod cluster;
mod combined;
mod decompose;
mod exact_solver;
mod external;
mod merge;
mod naive_upgma;
mod rt;
mod state;
mod utils;

use clap::{Parser, Subcommand};
use cluster::GCMStep;

use ordered_float::NotNan;

use std::path::PathBuf;
use tracing::info;

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
        #[clap(short, long, arg_enum, default_value_t = GCMStep::Auto)]
        tracer: GCMStep,
        #[clap(short, long, multiple_values = true)]
        weights: Vec<NotNan<f64>>,
        #[clap(short, long)]
        output: PathBuf,
    },

    Slice {
        #[clap(short, long)]
        input: PathBuf,
        #[clap(short, long)]
        tree: PathBuf,
        #[clap(short, long, value_parser = parse_axb, default_value = "10x200")]
        glues: (usize, usize),
        #[clap(short, long)]
        outdir: PathBuf,
        #[clap(short = 'c', long)]
        max_count: Option<usize>,
        #[clap(short = 's', long)]
        max_size: Option<usize>,
    },
}

fn parse_axb(s: &str) -> Result<(usize, usize), String> {
    let mut parts = s.split('x');
    let a = parts
        .next()
        .ok_or("missing first argument before 'x'".to_string())?;
    let b = parts
        .next()
        .ok_or("missing second argument before 'x'".to_string())?;
    let a_u = a
        .parse()
        .or_else(|_| Err("Cannot parse first argument into int"))?;
    let b_u = b
        .parse()
        .or_else(|_| Err("Cannot parse second argument into int"))?;
    Ok((a_u, b_u))
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
            tracer,
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
            combined::oneshot_merge_alignments(&input, &glues, tracer, &w, &output)
                .expect("Failed to merge alignments");
        }
        SubCommand::Slice {
            input,
            tree,
            glues,
            outdir,
            max_count,
            max_size,
        } => {
            info!("Analysis: slicing unaligned sequences.");
            combined::oneshot_slice_sequences(&input, &tree, glues, max_count, max_size, &outdir)?;
        }
    }
    info!("Total elapsed time: {:?}", now.elapsed());
    Ok(())
}
