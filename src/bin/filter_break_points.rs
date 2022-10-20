use std::path::PathBuf;

use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Alignments in the SAM file format.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Minimum number of breakpoints in [chunk] bp.
    #[arg(short, long, default_value_t = 10)]
    min_clip_reads: usize,
    /// Window to calculate the breakpoints.
    #[arg(short, long, default_value_t = 10)]
    window_size: usize,
}

const MARGIN: usize = 5;
use std::io::*;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let alignments = std::fs::File::open(&args.alignments).map(BufReader::new)?;
    let mut has_flush_myname = false;
    let mut sam_records = vec![];
    for line in alignments.lines().filter_map(|l| l.ok()) {
        if line.starts_with('@') {
            println!("{line}");
        } else {
            if !has_flush_myname {
                has_flush_myname = true;
                let filename = args.alignments.as_os_str().to_str().unwrap();
                let parameters = format!(
                    "--alignments {} --min_clip_reads {} --windows_size {}",
                    filename, args.min_clip_reads, args.window_size
                );
                println!(
                    "@PG\tID:filter_large_indel\tPN:mito_check\tCL:filter_break_points {parameters}"
                );
            }
            let sam = match bio_utils::sam::Sam::new(&line) {
                Some(res) => res,
                None => continue,
            };
            sam_records.push((line, sam));
        }
    }
    let breakpoints = survey_break_points(&sam_records, args.min_clip_reads, args.window_size);
    let range = args.window_size * MARGIN;
    sam_records.retain(|(_, sam)| {
        let (start, end) = sam.get_range();
        breakpoints.iter().any(|&position| {
            let start = start.max(position) - start.min(position);
            let end = end.max(position) - end.min(position);
            start.min(end) < range
        })
    });
    for (line, _) in sam_records.iter() {
        println!("{line}");
    }
    Ok(())
}

//Return a set of bp-points
fn survey_break_points(
    records: &[(String, bio_utils::sam::Sam)],
    min_size: usize,
    window_size: usize,
) -> Vec<usize> {
    let mut num_of_stst_reads = vec![];
    for (_, record) in records.iter() {
        let (start, stop) = record.get_range();
        let (start, stop) = (start / window_size, stop / window_size);
        if num_of_stst_reads.len() <= stop {
            let len = stop - num_of_stst_reads.len() + 1;
            num_of_stst_reads.extend(std::iter::repeat(0).take(len));
        }
        num_of_stst_reads[start] += 1;
        num_of_stst_reads[stop] += 1;
    }
    for (i, num) in num_of_stst_reads.iter().enumerate() {
        eprintln!("{}\t{num}", i * window_size);
    }
    num_of_stst_reads
        .iter()
        .enumerate()
        .filter_map(|(i, &num)| (min_size < num).then_some(i * window_size))
        .collect()
}
