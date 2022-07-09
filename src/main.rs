mod aln;
mod cluster;
mod combined;
mod exact_solver;
mod external;
mod merge;
mod naive_upgma;
mod progressive;
mod scorer;
mod state;
mod utils;

use clap::{Parser, Subcommand};
use cluster::GCMStep;
use ordered_float::NotNan;
use std::path::PathBuf;
use tracing::info;

use crate::{combined::autofix_input_constraints, scorer::oneshot_score_alignment};

#[derive(Parser, Debug, Hash, PartialEq)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
}

#[derive(Subcommand, Debug, PartialEq, Hash)]
enum SubCommand {
    /// Run GCM using existing subset and glue alignments
    Merge {
        /// Subset alignments
        #[clap(short, long, multiple_values = true)]
        input: Vec<PathBuf>,
        /// Glue alignments
        #[clap(short, long, multiple_values = true)]
        glues: Vec<PathBuf>,
        /// Tracing strategy
        #[clap(short, long, arg_enum, default_value_t = GCMStep::Auto)]
        tracer: GCMStep,
        /// Optional weights to the glues, same length as glue alignments
        #[clap(short, long, multiple_values = true)]
        weights: Vec<NotNan<f64>>,
        /// Output merged alignment path
        #[clap(short, long)]
        output: PathBuf,
    },

    /// Slice unaligned sequences into unaligned subsets and glues
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

    DebugImprove {
        /// Subset alignments
        #[clap(short, long, multiple_values = true)]
        input: Vec<PathBuf>,
        /// Serialized Graph
        #[clap(short, long)]
        graph: PathBuf,
        /// Serialized Trace
        #[clap(short, long)]
        trace: PathBuf,
        /// Output merged alignment path
        #[clap(short, long)]
        output: PathBuf,
    },

    DebugScore {
        #[clap(short, long, multiple_values = true)]
        input: Vec<PathBuf>,
        #[clap(short, long)]
        product: PathBuf,
        /// Serialized Graph
        #[clap(short, long)]
        graph: PathBuf,
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
        .map_err(|_| "Cannot parse first argument into int")?;
    let b_u = b
        .parse()
        .map_err(|_| "Cannot parse second argument into int")?;
    Ok((a_u, b_u))
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    use std::time::Instant;
    let now = Instant::now();
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .init();
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
        SubCommand::DebugImprove {
            mut input,
            graph,
            trace,
            output,
        } => {
            info!("Analysis: improving alignment.");
            autofix_input_constraints(&mut input);
            combined::oneshot_optimize_trace(&input, &graph, &trace, &output)?;
        }
        SubCommand::DebugScore {
            mut input,
            product,
            graph,
        } => {
            info!("Analysis: scoring alignment.");
            autofix_input_constraints(&mut input);
            let score = oneshot_score_alignment(&input, &graph, &product)?;
            println!("{}", serde_json::to_string_pretty(&score)?);
        }
    }
    info!("Total elapsed time: {:?}", now.elapsed());
    Ok(())
}
