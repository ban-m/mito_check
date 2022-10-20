use std::path::PathBuf;

use clap::Parser;

/// Enumerate break points in TSV format.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Alignments in Last's TSV format.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Minimum required match length.
    #[arg(short, long, default_value_t = 2000)]
    min_match_len: usize,
}

use std::io::*;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    println!("Contig\tStart\tEnd");
    let alignments = std::fs::File::open(&args.alignments)
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#'))
        .filter_map(|alignment| match_length_filter(&alignment, args.min_match_len));
    for (contig, start, end) in alignments {
        println!("{contig}\t{start}\t{end}");
    }
    Ok(())
}

fn match_length_filter(alignment: &str, min_length: usize) -> Option<(String, usize, usize)> {
    let alignment: Vec<_> = alignment.split('\t').collect();
    let match_len: usize = alignment[11]
        .split(',')
        .filter_map(|ops| ops.parse::<usize>().ok())
        .sum();
    let contig_name = alignment[6].to_string();
    let start: usize = alignment[7].parse().unwrap();
    let aln_size: usize = alignment[8].parse().unwrap();
    let strand = alignment[9] == "+";
    let seq_len: usize = alignment[10].parse().unwrap();
    let (start, end) = match strand {
        true => (start, start + aln_size),
        false => (seq_len - start - aln_size, seq_len - start),
    };
    if min_length < match_len && start != 0 && end != seq_len {
        Some((contig_name, start, end))
    } else {
        None
    }
}
