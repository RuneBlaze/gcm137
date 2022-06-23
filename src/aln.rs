use std::path::PathBuf;

use seq_io::fasta::{Reader, RefRecord};

pub trait AlnProcessor {
    type Output;
    fn on_record(&mut self, record: &RefRecord) -> anyhow::Result<()>;
    fn take(&mut self) -> Self::Output;
}

pub fn process_aln2<A: AlnProcessor + Default, B: AlnProcessor + Default>(
    infile: &PathBuf,
) -> anyhow::Result<(A::Output, B::Output)> {
    let mut p1 = A::default();
    let mut p2 = B::default();
    let mut reader = Reader::from_path(infile)?;
    while let Some(s) = reader.next() {
        let r = s?;
        p1.on_record(&r)?;
        p2.on_record(&r)?;
    }
    Ok((p1.take(), p2.take()))
}

pub fn process_aln1<A: AlnProcessor + Default>(infile: &PathBuf) -> anyhow::Result<A::Output> {
    let mut p = A::default();
    let mut reader = Reader::from_path(infile)?;
    while let Some(s) = reader.next() {
        let r = s?;
        p.on_record(&r)?;
    }
    Ok(p.take())
}
