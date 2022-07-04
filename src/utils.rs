use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};

use itertools::Itertools;
use rand::{prelude::ThreadRng, Rng};
use seq_io::{fasta::RefRecord, BaseRecord};

use crate::aln::AlnProcessor;

/// A streaming sequence sampler.
pub struct SequenceSampler {
    pub rng: ThreadRng,
    pub names: Vec<Vec<u8>>,
    pub records: Vec<Vec<u8>>,
    pub max_capacity: Option<usize>,
    pub i: usize,
}

fn degap(record: &RefRecord) -> Vec<u8> {
    let mut degapped = Vec::new();
    for b in record.seq().iter() {
        if *b == b'-' {
            continue;
        }
        degapped.push(*b);
    }
    degapped
}

impl SequenceSampler {
    pub fn new(max_capacitiy: Option<usize>) -> Self {
        Self {
            rng: rand::thread_rng(),
            names: Vec::new(),
            records: vec![],
            max_capacity: max_capacitiy,
            i: 0,
        }
    }

    pub fn is_full(&self) -> bool {
        self.max_capacity.map_or(false, |b| self.i >= b)
    }

    pub fn see(&mut self, record: &RefRecord) {
        match self.max_capacity {
            None => {
                self.names.push(record.head().to_vec());
                self.records.push(degap(record));
            }
            Some(s) => {
                if self.i < s {
                    self.names.push(record.head().to_vec());
                    self.records.push(degap(record));
                } else {
                    let j = self.rng.gen_range(0..(self.i + 1));
                    if j < s {
                        self.names[j] = record.head().to_vec();
                        self.records[j] = degap(record);
                    }
                }
            }
        }
        self.i += 1;
    }

    pub fn dump(&self, outfile: &PathBuf) -> anyhow::Result<()> {
        let mut writer = BufWriter::new(std::fs::File::create(outfile)?);
        for (name, seq) in self.names.iter().zip(self.records.iter()) {
            writer.write_all(b">")?;
            writer.write_all(name)?;
            writer.write_all(b"\n")?;
            seq.chunks(60).for_each(|chunk| {
                writer.write_all(chunk).unwrap();
                writer.write_all(b"\n").unwrap();
            });
        }
        Ok(())
    }
}

impl AlnProcessor for SequenceSampler {
    type Output = ();

    fn on_record(&mut self, record: &RefRecord) -> anyhow::Result<()> {
        self.see(record);
        Ok(())
    }

    fn take(&mut self) -> Self::Output {}
}

pub fn sample_seqs_in_place(
    src: &PathBuf,
    num_seqs: usize,
    out: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut sampler = SequenceSampler::new(Some(num_seqs));
    let mut reader = seq_io::fasta::Reader::from_path(src)?;
    while let Some(result) = reader.next() {
        let rec = result?;
        sampler.see(&rec);
    }
    let mut writer = BufWriter::new(std::fs::File::create(out)?);
    for (name, seq) in sampler.names.iter().zip(sampler.records.iter()) {
        writer.write_all(b">")?;
        writer.write_all(name)?;
        writer.write_all(b"\n")?;
        seq.chunks(60).for_each(|chunk| {
            writer.write_all(chunk).unwrap();
            writer.write_all(b"\n").unwrap();
        });
    }
    Ok(())
}
