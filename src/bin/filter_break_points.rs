use std::path::PathBuf;

use clap::Parser;

/// Filter sam file to detect linear structures.
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
    /// Merge window.
    #[arg(short, long, default_value_t = 10000)]
    merge_window: usize,
    /// Squish window
    #[arg(short, long, default_value_t = 2000)]
    squish_window: usize,
}

const MARGIN: usize = 5;
use std::io::*;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    if args.alignments.as_os_str() == "-" {
        let stdin = std::io::stdin();
        let alignments = BufReader::new(stdin.lock()).lines().filter_map(|l| l.ok());
        flush_break_reads(alignments, &args);
    } else {
        let alignments = std::fs::File::open(&args.alignments).map(BufReader::new)?;
        let alignments = alignments.lines().filter_map(|l| l.ok());
        flush_break_reads(alignments, &args);
    }
    Ok(())
}

fn flush_break_reads<I: std::iter::Iterator<Item = String>>(alignments: I, args: &Args) {
    let mut has_flush_myname = false;
    let mut sam_records = vec![];
    for line in alignments {
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
    let mut breakpoints = survey_break_points(&sam_records, args.min_clip_reads, args.window_size);
    breakpoints
        .values_mut()
        .for_each(|positions| *positions = squish_near_breaks(positions, args.squish_window));
    let breakpoints = split_breakpoints(&breakpoints, args.merge_window);
    for (id, (iso, pair)) in breakpoints.iter() {
        for pos in iso.iter() {
            eprintln!("{id}\tIso\t{pos}");
        }
        for (start, end) in pair.iter() {
            eprintln!("{id}\tPair\t{start}\t{end}");
        }
    }
    let range = args.window_size * MARGIN;
    sam_records.retain(|(_, sam)| is_in_breakpoints(&breakpoints, sam, range));
    for (line, _) in sam_records.iter() {
        println!("{line}");
    }
}

fn squish_near_breaks(positions: &[usize], window_size: usize) -> Vec<usize> {
    if positions.is_empty() {
        return Vec::new();
    }
    let mut squished = vec![];
    let mut current = positions[0];
    let mut sums_so_far = current;
    let mut len = 1;
    for &pos in positions.iter().skip(1) {
        if current + window_size < pos {
            squished.push(sums_so_far / len);
            sums_so_far = pos;
            len = 1;
        } else {
            sums_so_far += pos;
            len += 1;
        }
        current = pos;
    }
    squished.push(sums_so_far / len);
    squished.pop();
    squished
}

type BreakPoints = HashMap<String, Breaks>;
type Breaks = (Vec<usize>, Vec<(usize, usize)>);
fn split_breakpoints(
    breakpoints: &HashMap<String, Vec<usize>>,
    merge_window: usize,
) -> BreakPoints {
    breakpoints
        .iter()
        .map(|(id, positions)| {
            let (mut isolated, mut paired) = (vec![], vec![]);
            let mut is_paired = false;
            for (i, &pos) in positions.iter().enumerate() {
                match positions.get(i + 1) {
                    Some(&next) if next < pos + merge_window => {
                        paired.push((pos, next));
                        is_paired = true;
                    }
                    Some(_) if is_paired => {
                        is_paired = false;
                    }
                    Some(_) => {
                        is_paired = false;
                        isolated.push(pos);
                    }
                    None if is_paired => {}
                    None => isolated.push(pos),
                }
            }
            (id.clone(), (isolated, paired))
        })
        .collect()
}

fn is_in_breakpoints(breakpoints: &BreakPoints, sam: &bio_utils::sam::Sam, range: usize) -> bool {
    let (isolated, paired) = breakpoints.get(sam.ref_name()).unwrap();
    let (start, end) = sam.get_range();
    let touch_isolated = isolated
        .iter()
        .any(|&pos| pos.max(start) - pos.min(start) < range || pos.max(end) - pos.min(end) < range);
    let contained_in_pair = paired
        .iter()
        .any(|&(l, r)| l.saturating_sub(range) < start && end < r + range);
    touch_isolated || contained_in_pair
}

use std::collections::HashMap;
//Return a set of bp-points
fn survey_break_points(
    records: &[(String, bio_utils::sam::Sam)],
    min_size: usize,
    window_size: usize,
) -> HashMap<String, Vec<usize>> {
    let mut num_of_stst_reads: HashMap<_, HashMap<_, usize>> = HashMap::new();
    for (_, record) in records.iter() {
        let (start, stop) = record.get_range();
        let (start, stop) = (start / window_size, stop / window_size);
        let ref_name = record.ref_name().to_string();
        let slot = num_of_stst_reads.entry(ref_name).or_default();
        *slot.entry(start).or_default() += 1;
        *slot.entry(stop).or_default() += 1;
    }
    num_of_stst_reads
        .into_iter()
        .map(|(id, counts)| {
            let mut positions: Vec<usize> = counts
                .into_iter()
                .filter_map(|(pos, count)| (min_size < count).then_some(window_size * pos))
                .collect();
            positions.sort_unstable();
            positions.retain(|&x| x != 0);
            (id, positions)
        })
        .collect()
}
