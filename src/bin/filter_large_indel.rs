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
    let mut has_flush_myname = false;
    for line in alignments.lines().filter_map(|l| l.ok()) {
        if line.starts_with('@') {
            println!("{line}");
        } else {
            if !has_flush_myname {
                has_flush_myname = true;
                let filename = args.alignments.as_os_str().to_str().unwrap();
                println!("@PG\tID:filter_large_indel\tPN:mito_check\tCL:filter_large_indel --alignments {} --min_sv_size {}", filename, args.min_sv_size);
            }
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
