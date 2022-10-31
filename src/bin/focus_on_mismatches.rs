use clap::Parser;
use std::path::PathBuf;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The reference sequence.
    #[arg(short, long)]
    reference: PathBuf,
    /// The assembled contigs.
    #[arg(short, long)]
    contigs: PathBuf,
    /// Alignments between the reference and the contigs.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Minimum crop range.
    #[arg(long, default_value_t = 20)]
    min_crop_len: usize,
    /// Maximum crop range.
    #[arg(long, default_value_t = 500)]
    max_crop_len: usize,
    #[arg(long)]
    prefix: PathBuf,
    /// Remove mismatches [Length]-bp near the start/end of the contig.
    #[arg(long, default_value_t = 1000)]
    safe_margin: usize,
}

use std::collections::HashMap;
use std::io::*;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let reference = bio_utils::fasta::parse_into_vec(&args.reference)?;
    let contigs = bio_utils::fasta::parse_into_vec(&args.contigs)?;
    let mut alignments: HashMap<_, Vec<_>> = HashMap::new();
    for aln in std::fs::File::open(&args.alignments)
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .filter_map(|l| bio_utils::paf::PAF::new(&l))
    {
        let query = aln.qname.clone();
        alignments.entry(query).or_default().push(aln);
    }
    let crop_regions = find_crop_regions(&alignments, &args);
    for crop_region in crop_regions.iter() {
        let (sid, (start, end)) = crop_region;
        let mut output = args.prefix.clone();
        output.push(format!("{sid}_{start}_{end}.fa"));
        let mut wtr = std::fs::File::create(output).map(BufWriter::new)?;
        let mut sequences: Vec<_> = reference
            .iter()
            .find_map(|refr| {
                if refr.id() == sid {
                    let seq = String::from_utf8_lossy(&refr.seq()[*start..*end]).to_string();
                    let id = format!("{}_{start}_{end}", refr.id());
                    Some(vec![(id, seq)])
                } else {
                    None
                }
            })
            .unwrap();
        for contig in contigs.iter() {
            let alns = match alignments.get(contig.id()) {
                Some(res) => res,
                None => continue,
            };
            sequences.extend(get_crop_region(contig, alns, crop_region));
        }
        for (id, seq) in sequences {
            if !seq.is_empty() {
                writeln!(&mut wtr, ">{id}\n{seq}")?;
            }
        }
    }
    Ok(())
}

use bio_utils::paf::PAF;
type Region = (String, (usize, usize));
fn find_crop_regions(alignments: &HashMap<String, Vec<PAF>>, arg: &Args) -> Vec<Region> {
    let position_of_variants = enumerate_mismatches(alignments, arg);
    position_of_variants
        .iter()
        .flat_map(|(id, poss)| {
            let regions = aggregate(poss, arg);
            regions.into_iter().map(|range| (id.clone(), range))
        })
        .collect()
}

fn aggregate(poss: &[usize], args: &Args) -> Vec<(usize, usize)> {
    if poss.is_empty() {
        return Vec::new();
    }
    let mut start_pos = poss[0];
    let mut current_pos = poss[0];
    let mut ranges = vec![];
    for &pos in poss.iter().skip(1) {
        if args.min_crop_len < pos - current_pos || args.max_crop_len < current_pos - start_pos {
            let range = (
                start_pos.saturating_sub(args.min_crop_len),
                current_pos + args.min_crop_len,
            );
            ranges.push(range);
            start_pos = pos;
            current_pos = pos;
        } else {
            current_pos = pos;
        }
    }
    ranges
}

fn enumerate_mismatches(
    alignments: &HashMap<String, Vec<PAF>>,
    arg: &Args,
) -> HashMap<String, Vec<usize>> {
    let mut found_mismatches: HashMap<_, Vec<_>> = HashMap::new();
    for aln in alignments.values().flat_map(|alns| alns.iter()) {
        let mut rpos = aln.tstart;
        let slot = found_mismatches.entry(aln.tname.clone()).or_default();
        let cigar = bio_utils::sam::parse_cigar_string(aln.get_tag("cg").unwrap().1);
        for op in cigar {
            use bio_utils::sam::Op;
            match op {
                Op::Match(l) | Op::Align(l) => {
                    rpos += l;
                }
                Op::Insertion(_) => {}
                Op::Deletion(l) => rpos += l,
                Op::Mismatch(l) => {
                    slot.push(rpos);
                    rpos += l;
                }
                _ => {}
            }
        }
    }
    let contig_lengths: HashMap<_, _> = alignments
        .values()
        .flat_map(|alns| alns.iter())
        .map(|aln| (aln.tname.to_string(), aln.tlen))
        .collect();
    for (contig_name, slot) in found_mismatches.iter_mut() {
        slot.sort_unstable();
        slot.dedup();
        let length = contig_lengths[contig_name];
        slot.retain(|&pos| arg.safe_margin < pos && pos < length.saturating_sub(arg.safe_margin));
    }
    found_mismatches
}

use bio_utils::fasta::Record;
fn get_crop_region(
    contig: &Record,
    alns: &[PAF],
    (ref_id, range): &Region,
) -> Vec<(String, String)> {
    let seq = contig.seq();
    let &(start, end) = range;
    let id = contig.id();
    alns.iter()
        .filter_map(|aln| {
            let out_of_range = aln.tend < start || end < aln.tstart;
            if &aln.tname != ref_id || out_of_range {
                None
            } else {
                let (qstart, qend, direction) = crop(aln, start, end);
                let seq = match direction {
                    true => String::from_utf8_lossy(&seq[qstart..qend]).to_string(),
                    false => String::from_utf8(bio_utils::revcmp(&seq[qstart..qend])).unwrap(),
                };
                let id = format!("{id}_{qstart}_{qend}_{direction}");
                Some((id, seq))
            }
        })
        .collect()
}

fn crop(alignment: &PAF, start: usize, end: usize) -> (usize, usize, bool) {
    let mut crop_start = None;
    let mut crop_end = None;
    let (mut rpos, mut qpos) = (alignment.tstart, 0);
    let cigar = bio_utils::sam::parse_cigar_string(alignment.get_tag("cg").unwrap().1);
    let mut ops = cigar.iter().flat_map(|op| {
        use bio_utils::sam::Op;
        let (op, len) = match op {
            Op::Align(l) | Op::Mismatch(l) | Op::Match(l) => (b'M', *l),
            Op::Insertion(l) | Op::SoftClip(l) | Op::HardClip(l) => (b'I', *l),
            Op::Deletion(l) => (b'D', *l),
            _ => panic!(),
        };
        std::iter::repeat(op).take(len)
    });
    while rpos < start {
        match ops.next() {
            Some(b'M') => {
                rpos += 1;
                qpos += 1
            }
            Some(b'D') => rpos += 1,
            Some(b'I') => qpos += 1,
            None => break,
            _ => panic!(),
        }
    }
    if rpos == start {
        crop_start = Some(qpos);
    }
    while rpos < end {
        match ops.next() {
            Some(b'M') => {
                rpos += 1;
                qpos += 1
            }
            Some(b'D') => rpos += 1,
            Some(b'I') => qpos += 1,
            None => break,
            _ => panic!(),
        }
    }
    if rpos == end {
        crop_end = Some(qpos);
    }
    // eprintln!(
    //     "{}\t{}\t{crop_start:?}\t{crop_end:?}",
    //     alignment.qname, alignment.relstrand
    // );
    if alignment.relstrand {
        let c_start = match crop_start {
            Some(res) => alignment.qstart + res,
            None => alignment.qstart,
        };
        let c_end = match crop_end {
            Some(res) => alignment.qstart + res,
            None => alignment.qend,
        };
        (c_start, c_end, true)
    } else {
        let c_start = match crop_end {
            Some(res) => alignment.qend - res,
            None => alignment.qstart,
        };
        let c_end = match crop_start {
            Some(res) => alignment.qend - res,
            None => alignment.qend,
        };
        (c_start, c_end, false)
    }
}
