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
    println!("Contig\tCStart\tCEnd\tRefr\tRStart\tREnd");
    let alignments = std::fs::File::open(&args.alignments)
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#'))
        .filter_map(|alignment| match_length_filter(&alignment, args.min_match_len));
    for ((contig, c_start, c_end), (refr, r_start, r_end)) in alignments {
        println!("{contig}\t{c_start}\t{c_end}\t{refr}\t{r_start}\t{r_end}");
    }
    Ok(())
}

type BreakPoint = (String, usize, usize);
fn match_length_filter(alignment: &str, min_length: usize) -> Option<(BreakPoint, BreakPoint)> {
    let alignment: Vec<_> = alignment.split('\t').collect();
    let match_len: usize = alignment[11]
        .split(',')
        .filter_map(|ops| ops.parse::<usize>().ok())
        .sum();
    let (ref_break_point, ref_len) = {
        let ref_name = alignment[1].to_string();
        let start: usize = alignment[2].parse().unwrap();
        let aln_size: usize = alignment[3].parse().unwrap();
        let is_forward = alignment[4] == "+";
        let size: usize = alignment[5].parse().unwrap();
        let break_point = match is_forward {
            true => (ref_name, start, start + aln_size),
            false => (ref_name, size - start - aln_size, size - start),
        };
        (break_point, size)
    };
    // Remove SV near the start/end of a contig.
    let is_refr_boundary = ref_break_point.1 < 100 || ref_len - ref_break_point.2 < 100;
    let (contig_break_point, seq_len) = {
        let contig_name = alignment[6].to_string();
        let start: usize = alignment[7].parse().unwrap();
        let aln_size: usize = alignment[8].parse().unwrap();
        let strand = alignment[9] == "+";
        let seq_len: usize = alignment[10].parse().unwrap();
        match strand {
            true => ((contig_name, start, start + aln_size), seq_len),
            false => (
                (contig_name, seq_len - start - aln_size, seq_len - start),
                seq_len,
            ),
        }
    };
    let is_contig_boundary = contig_break_point.1 < 100 || seq_len - contig_break_point.2 < 100;
    if min_length < match_len && !is_contig_boundary && !is_refr_boundary {
        Some((contig_break_point, ref_break_point))
    } else {
        None
    }
}
