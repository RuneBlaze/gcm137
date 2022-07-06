use crate::{cluster::{GCMStep, ClusteringResult}, exact_solver::sw_algorithm, progressive::iterative_refinement, merge::load_graph};
use anyhow::Ok;
use ogcat::ogtree::{self, TreeCollection};
use ordered_float::NotNan;
use rand::prelude::SliceRandom;
use seq_io::{
    fasta::{OwnedRecord, Reader},
    BaseRecord,
};
use std::{
    fs::create_dir_all,
    io::{BufWriter, Write},
    path::PathBuf,
};
use tracing::{debug, warn, info};

use crate::{
    aln::AlnProcessor,
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
    tracer_mode: GCMStep,
    weights: &Option<Vec<NotNan<f64>>>,
    outpath: &PathBuf,
) -> anyhow::Result<()> {
    let mut state = state_from_constraints(constraints)?;
    debug!("Constructed state from constraints");
    let graph = build_graph(&mut state, glues, weights).unwrap();
    debug!("Built alignment graph.");
    let mut res = if constraints.len() == 2 && tracer_mode != GCMStep::Upgma {
        debug!("Running Smith-Waterman, solving MWT-AM exactly.");
        sw_algorithm(&graph, &state)
    } else {
        debug!("Running UPGMA heuristic for solving MWT-AM.");
        naive_upgma(&graph, &state)
    };
    debug!("Clustered/Traced alignment graph.");
    iterative_refinement(&state, &graph, &mut res);
    debug!("Finished iterative refinement.");
    let frames = build_frames(&state, &res);
    debug!("Flushing merged alignments...");
    merge_alignments_from_frames(constraints, &frames, outpath)?;
    Ok(())
}

pub fn oneshot_stitch_alignments(
    constraints: &[PathBuf],
    _outpath: &PathBuf,
) -> anyhow::Result<()> {
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
    todo!(); // honestly I don't know if I should implement this or not
    Ok(())
}

pub fn oneshot_optimize_trace(
    constraints: &[PathBuf],
    graph_path: &PathBuf,
    trace_path: &PathBuf,
    outpath: &PathBuf,
) -> anyhow::Result<()> {
    let mut state = state_from_constraints(constraints)?;
    debug!("Constructed state from constraints");
    let mut graph = load_graph(&state, graph_path)?;
    let mut trace = ClusteringResult::from_plaintext(trace_path, &graph)?;
    let before_score = trace.mwt_am_score(&state, &graph);
    iterative_refinement(&state, &graph, &mut trace);
    let after_score = trace.mwt_am_score(&state, &graph);
    info!("Optimized trace: {:.2} -> {:.2}, ({:.2}% increase)", before_score, after_score, (after_score - before_score) / before_score * 100.0);
    let frames = build_frames(&state, &trace);
    debug!("Flushing merged alignments...");
    merge_alignments_from_frames(constraints, &frames, outpath)?;
    Ok(())
}

pub fn oneshot_slice_sequences(
    input: &PathBuf,
    tree: &PathBuf,
    glues: (usize, usize),
    max_count: Option<usize>,
    max_size: Option<usize>,
    outdir: &PathBuf,
) -> anyhow::Result<()> {
    let mut rng = rand::thread_rng();
    let collection = TreeCollection::from_newick(tree).expect("Failed to read tree");
    let decomp = ogtree::centroid_edge_decomp(&collection.trees[0], &max_count, &max_size);
    let labels = ogtree::cuts_to_subsets(&collection.trees[0], &decomp);
    let ts = &collection.taxon_set;
    let mut reader = Reader::from_path(input)?;
    let mut subsets: Vec<Vec<OwnedRecord>> = vec![Vec::new(); decomp.len()];
    while let Some(s) = reader.next() {
        let r = s?;
        let head = String::from_utf8_lossy(r.head());
        let id = ts.retrieve(&head);
        subsets[labels[id]].push(r.to_owned_record());
    }
    create_dir_all(outdir)?;
    let (glue_num, glue_size) = glues;
    let mut constraints_path = outdir.clone();
    constraints_path.push("constraints");
    create_dir_all(&constraints_path)?;
    let mut glues_path = outdir.clone();
    glues_path.push("glues");
    create_dir_all(&glues_path)?; // oh my god this is so ugly
    for (i, c) in subsets.iter().enumerate() {
        let mut cp = constraints_path.clone();
        cp.push(format!("constraint_{}.unaln.fa", i));
        let file = std::fs::File::create(cp)?;
        let mut writer = BufWriter::new(file);
        for r in c {
            writer.write_all(b">")?;
            writer.write_all(r.head())?;
            writer.write_all(b"\n")?;
            r.seq.chunks(60).try_for_each(|chunk| {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
                Ok(())
            })?;
        }
    }
    let k = decomp.len();
    let mut sample_size = glue_size / k;
    let max_seqset_size = subsets.iter().map(|it| it.len()).max().unwrap();
    if sample_size > max_seqset_size {
        warn!(
            "Sample size is larger than largest subset size. Setting sample size to {}",
            max_seqset_size
        );
        sample_size = max_seqset_size;
    }
    for i in 0..glue_num {
        let mut gp = glues_path.clone();
        gp.push(format!("glue_{}.unaln.fa", i));
        let file = std::fs::File::create(gp)?;
        let mut writer = BufWriter::new(file);
        for c in &subsets {
            c.choose_multiple(&mut rng, sample_size).try_for_each(|r| {
                writer.write_all(b">")?;
                writer.write_all(r.head())?;
                writer.write_all(b"\n")?;
                r.seq.chunks(60).try_for_each(|chunk| {
                    writer.write_all(chunk)?;
                    writer.write_all(b"\n")?;
                    Ok(())
                })?;
                Ok(())
            })?;
        }
    }
    Ok(())
}
