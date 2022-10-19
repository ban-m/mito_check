use std::path::PathBuf;

use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Alignments in the SAM file format.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Minimum length of the insertion/deletion to be retained.
    #[arg(short, long, default_value_t = 500)]
    min_sv_size: usize,
}

use std::io::*;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let alignments = std::fs::File::open(&args.alignments).map(BufReader::new)?;
    for line in alignments.lines().filter_map(|l| l.ok()) {
        if line.starts_with('@') {
            println!("{line}");
        } else {
            let sam = match bio_utils::sam::Sam::new(&line) {
                Some(res) => res,
                None => continue,
            };
            use bio_utils::sam::Op;
            let cigar = sam.cigar();
            let max_indel = cigar
                .iter()
                .map(|op| match op {
                    Op::Deletion(l) | Op::Insertion(l) => *l,
                    _ => 0,
                })
                .max()
                .unwrap_or(0);
            if args.min_sv_size < max_indel {
                println!("{line}");
            }
        }
    }
    Ok(())
}
