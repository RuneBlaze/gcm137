use std::path::PathBuf;
use tokio::{
    fs::{rename, File},
    io::BufWriter,
    process::Command,
};

/// runs mafft-linsi
async fn request_alignment(
    in_path: &PathBuf,
    out_path: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
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
        .spawn()?;
    let mut temp_out = out_path.clone();
    temp_out.set_extension("temp");
    let mut stdout = child.stdout.take().unwrap();
    let f = File::create(&temp_out).await?;
    let mut writer = BufWriter::new(f);
    tokio::io::copy(&mut stdout, &mut writer).await?;
    rename(&temp_out, &out_path).await?;
    Ok(())
}
