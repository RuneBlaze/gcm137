

use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use std::path::PathBuf;
use std::{
    error::Error,
};

use crate::{
    merge::{build_frames, build_graph, merge_alignments_from_frames, state_from_constraints},
    naive_upgma::naive_upgma,
};

pub fn oneshot_merge_alignments(
    constraints: &[PathBuf],
    glues: &[PathBuf],
    outpath: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut state = state_from_constraints(constraints)?;
    println!("Constructed state from constraints");
    let graph = build_graph(&mut state, glues).unwrap();
    println!("Built graph");
    let res = naive_upgma(&graph, &state);
    println!("Built clusters: {}", res.clusters.len());
    // res.check_validity();
    // println!("Checked validity");
    // println!("clusters {:?}", res.clusters);
    // assert_ne!(1,1);
    let frames = build_frames(&state, &res);
    let _frame_lengths = frames
        .iter()
        .map(|f| {
            return f.iter().sum::<u32>() + f.len() as u32;
        })
        .collect_vec();
    // println!("built frames : {:?} {:?}, {:?}",
    // frame_lengths, frames.iter().map(|it| it.len()).collect_vec(), &state.column_counts);
    merge_alignments_from_frames(constraints, &frames, outpath)?;
    Ok(())
}
