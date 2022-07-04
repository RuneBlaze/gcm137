use std::fs::File;
use std::os::unix::prelude::{FromRawFd, IntoRawFd};
use std::process::Stdio;
use std::{fs::rename, path::PathBuf, process::Command};

pub fn request_alignment(in_path: &PathBuf, out_path: &PathBuf) -> anyhow::Result<()> {
    let mut temp_out = out_path.clone();
    temp_out.set_extension("temp");
    let fd = File::create(&temp_out)?.into_raw_fd();
    let mut child = Command::new("mafft")
        .arg("--localpair")
        .arg("--maxiterate")
        .arg("1000")
        .arg("--ep")
        .arg("0.123")
        .arg("--quiet")
        .arg("--thread")
        .arg("2")
        .arg(in_path)
        .stdout(unsafe { Stdio::from_raw_fd(fd) })
        .spawn()?;
    child.wait()?;
    rename(&temp_out, &out_path)?;
    Ok(())
}
