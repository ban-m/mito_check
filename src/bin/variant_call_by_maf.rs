use std::path::PathBuf;

use clap::Parser;

/// Call variants by MAF.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Alignments in last MAF format.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Minimum required match length.
    #[arg(short, long, default_value_t = 2000)]
    min_match_len: usize,
}

fn main() -> std::io::Result<()> {
    let _ = Args::parse();
    Ok(())
}
