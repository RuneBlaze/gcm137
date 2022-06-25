use crate::{cluster::ClusteringResult, exact_solver::sw_algorithm};
use anyhow::Ok;
use itertools::Itertools;
use ordered_float::NotNan;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use seq_io::fasta::Reader;
use tracing::{debug, error, info, span, warn, Level};

use std::error::Error;
use std::path::PathBuf;

use crate::{
    aln::AlnProcessor,
    exact_solver,
    merge::{
        build_frames, build_graph, merge_alignments_from_frames, state_from_constraints,
        StateFromConstraints,
    },
    naive_upgma::naive_upgma,
    utils::SequenceSampler,
};

#[tracing::instrument]
pub fn oneshot_merge_alignments(
    constraints: &[PathBuf],
    glues: &[PathBuf],
    weights: &Option<Vec<NotNan<f64>>>,
    outpath: &PathBuf,
) -> anyhow::Result<()> {
    let mut state = state_from_constraints(constraints)?;
    debug!("Constructed state from constraints");
    let graph = build_graph(&mut state, glues, weights).unwrap();
    debug!("Built alignment graph.");
    let res = if constraints.len() == 2 {
        debug!("Running SW algorithm, solving MWT-AM exactly.");
        sw_algorithm(&graph, &state)
    } else {
        debug!("Running UPGMA heuristic for solving MWT-AM.");
        naive_upgma(&graph, &state)
    };
    debug!("Clustered/Traced alignment graph.");
    let frames = build_frames(&state, &res);
    // println!("Built frames. {:?}", frames[0]);
    debug!("Flushing merged alignments...");
    merge_alignments_from_frames(constraints, &frames, outpath).unwrap();
    Ok(())
}

pub fn oneshot_stitch_alignments(constraints: &[PathBuf], outpath: &PathBuf) -> anyhow::Result<()> {
    let mut p = StateFromConstraints::default();
    for aln in constraints {
        // cid : constraint id
        let mut s = SequenceSampler::new(Some(150 / constraints.len()));
        let mut reader = Reader::from_path(aln)?;
        while let Some(result) = reader.next() {
            let rec = result?;
            p.on_record(&rec)?;
            s.on_record(&rec)?;
        }
        p.next_aln();
    }
    todo!();
    Ok(())
}
