use anyhow::Ok;
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use seq_io::fasta::Reader;

use std::error::Error;
use std::path::PathBuf;

use crate::{
    aln::AlnProcessor,
    merge::{
        build_frames, build_graph, merge_alignments_from_frames, state_from_constraints,
        StateFromConstraints,
    },
    naive_upgma::naive_upgma,
    utils::SequenceSampler,
};

pub fn oneshot_merge_alignments(
    constraints: &[PathBuf],
    glues: &[PathBuf],
    outpath: &PathBuf,
) -> anyhow::Result<()> {
    let mut state = state_from_constraints(constraints)?;
    println!("Constructed state from constraints");
    let graph = build_graph(&mut state, glues).unwrap();
    println!("Built graph");
    let res = naive_upgma(&graph, &state);
    println!("Built clusters: {}", res.clusters.len());
    let frames = build_frames(&state, &res);
    let _frame_lengths = frames
        .iter()
        .map(|f| {
            return f.iter().sum::<u32>() + f.len() as u32;
        })
        .collect_vec();
    merge_alignments_from_frames(constraints, &frames, outpath).unwrap();
    Ok(())
}

pub fn oneshot_stitch_alignments(constraints: &[PathBuf], outpath: &PathBuf) -> anyhow::Result<()> {
    let mut p = StateFromConstraints::default();
    let mut s = SequenceSampler::new(Some(150 / constraints.len()));
    for aln in constraints {
        // cid : constraint id
        let mut reader = Reader::from_path(aln)?;
        while let Some(result) = reader.next() {
            let rec = result?;
            p.on_record(&rec)?;
            s.on_record(&rec)?;
        }
        p.next_aln();
    }
    Ok(())
}
