use std::path::PathBuf;

use seq_io::fasta::RefRecord;

pub trait AlnProcessor {
    fn on_record(&mut self, record : &RefRecord);
}

pub fn process_aln<A : AlnProcessor + Default, B : AlnProcessor + Default>(infile : &PathBuf) -> (A, B) {
    let mut p1 = A::default();
    let mut p2 = B::default();
    todo!();
}