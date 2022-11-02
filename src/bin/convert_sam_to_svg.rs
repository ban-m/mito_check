use clap::Parser;
use std::path::PathBuf;

/// Produce SVG file from SAM file.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Alignments in the SAM file format. If -, use stdin.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Contigs.
    #[arg(short, long)]
    contigs: PathBuf,
    /// Output path
    #[arg(short, long)]
    output: PathBuf,
    // To re-calculate the alignment again.
    // #[arg(short, long)]
    // re_align: bool,
    /// To squish matches/deletions smaller than LEN.
    #[arg(short, long, default_value_t = 7)]
    squish: usize,
    /// Target region ex. chr1:10000-1500000.
    #[arg(short, long)]
    target: Option<String>,
}

use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let contigs = bio_utils::fasta::parse_into_vec(&args.contigs)?;
    let selection = match &args.target {
        Some(res) => Selection::new(res),
        None => Selection::from_contigs(&contigs),
    };
    let (_, mut samfile) = if args.alignments.as_os_str() == "-" {
        let stdin = std::io::stdin();
        let stdin = BufReader::new(stdin.lock()).lines().filter_map(|l| l.ok());
        parse_sam_file(stdin)
    } else {
        let file = std::fs::File::open(&args.alignments).map(BufReader::new)?;
        let file = file.lines().filter_map(|l| l.ok());
        parse_sam_file(file)
    };
    samfile.retain(|record| selection.is_in(record.r_name(), record.get_range()));
    let contigs: HashMap<String, _> = contigs
        .into_iter()
        .filter(|r| selection.has_contig(r.id()))
        .map(|r| (r.id().to_string(), r))
        .collect();
    let tidy_alignments = convert_tidy(samfile, &contigs, &args);
    flush_alignments(&tidy_alignments, &selection, &args);
    Ok(())
}

#[derive(Debug, Clone)]
struct Selection {
    selects: Vec<(String, usize, usize)>,
}

impl Selection {
    fn has_contig(&self, id: &str) -> bool {
        self.selects.iter().any(|x| x.0 == id)
    }
    fn new(select: &str) -> Self {
        let (name, range) = select.split_once(':').unwrap();
        let (start, end) = range.split_once('-').unwrap();
        let start: usize = start.parse().unwrap();
        let end: usize = end.parse().unwrap();
        let select = (name.to_string(), start, end);
        Self {
            selects: vec![select],
        }
    }
    fn from_contigs(contigs: &[bio_utils::fasta::Record]) -> Self {
        let selects: Vec<_> = contigs
            .iter()
            .map(|contig| (contig.id().to_string(), 0, contig.seq().len()))
            .collect();
        Self { selects }
    }
    fn is_in(&self, id: &str, (start, end): (usize, usize)) -> bool {
        self.selects
            .iter()
            .any(|&(ref ref_id, s, e)| id == ref_id && s <= start && end <= e)
    }
}

#[derive(Debug, Clone)]
struct TidyAlignments {
    // base pair.
    start: usize,
    end: usize,
    // Plot row.
    row: usize,
    operations: Vec<(usize, char)>,
    // Head/Tail clippings.
    clips: (usize, usize),
    contig: String,
    read_id: String,
}

fn convert_tidy(
    mut alignments: Vec<sam::Sam>,
    contigs: &HashMap<String, bio_utils::fasta::Record>,
    args: &Args,
) -> Vec<(String, Vec<TidyAlignments>)> {
    let mut tidyalignments: HashMap<String, Vec<_>> =
        contigs.keys().cloned().map(|k| (k, vec![])).collect();
    // Sort the alignments by the start position.
    alignments.sort_by_cached_key(|aln| {
        let range = aln.get_range();
        (aln.ref_name().to_string(), range)
    });
    // Contig ID -> Row Number -> Rightmost position.
    let mut row_fronteir: HashMap<String, Vec<usize>> =
        contigs.keys().cloned().map(|k| (k, vec![])).collect();
    for aln in alignments.iter() {
        let mut tidied = to_tidy_alignments(aln, contigs, args);
        let fronteir = row_fronteir.get_mut(&tidied.contig).unwrap();
        let row = fronteir.iter().position(|&n| n <= tidied.start);
        if let Some(row) = row {
            tidied.row = row;
            fronteir[row] = tidied.end;
        } else {
            tidied.row = fronteir.len();
            fronteir.push(tidied.end);
        }
        tidyalignments.get_mut(&tidied.contig).unwrap().push(tidied);
    }
    tidyalignments.into_iter().collect()
}

use bio_utils::sam;
fn to_tidy_alignments(
    aln: &bio_utils::sam::Sam,
    _contigs: &HashMap<String, bio_utils::fasta::Record>,
    args: &Args,
) -> TidyAlignments {
    let (start, end) = aln.get_range();
    let mut cigar = aln.cigar();
    let tail_clip = match cigar.last().cloned() {
        Some(sam::Op::HardClip(l)) => {
            cigar.pop();
            l
        }
        _ => 0,
    };
    let head_clip = match cigar.first().cloned() {
        Some(sam::Op::HardClip(l)) => {
            cigar.reverse();
            cigar.pop();
            cigar.reverse();
            l
        }
        _ => 0,
    };
    let clips = (head_clip, tail_clip);
    // if args.re_align {
    //     let contig = contigs.get(aln.ref_name()).unwrap();
    //     cigar = re_align(aln, contig, &cigar);
    // }
    let ops = to_op_runs(&cigar);
    let ops = compress_small_matches(&ops, args.squish);
    let ops = compress_small_deletions(&ops, args.squish);
    TidyAlignments {
        start,
        end,
        row: 0,
        operations: ops,
        clips,
        contig: aln.r_name().to_string(),
        read_id: aln.q_name().to_string(),
    }
}

// const RADIUS: usize = 200;
// const PARAMS: (i32, i32, i32, i32) = (5, -30, -35, -1);

// fn re_align(
//     aln: &bio_utils::sam::Sam,
//     contig: &bio_utils::fasta::Record,
//     cigar: &[sam::Op],
// ) -> Vec<sam::Op> {
//     let seq = aln.seq();
//     let len: usize = cigar
//         .iter()
//         .map(|op| match op {
//             sam::Op::Match(l)
//             | sam::Op::Mismatch(l)
//             | sam::Op::Align(l)
//             | sam::Op::Insertion(l) => *l,
//             _ => 0,
//         })
//         .sum();
//     assert_eq!(seq.len(), len);
//     let (start, end) = aln.get_range();
//     let contig = &contig.seq()[start..end];
//     let kiley_op = to_kiley_ops(cigar);
//     let kiley_op = if aln.is_forward() {
//         let seq = seq.as_bytes();
//         kiley::bialignment::guided::global_guided(contig, seq, &kiley_op, RADIUS, PARAMS).1
//     } else {
//         let seq = bio_utils::revcmp(seq.as_bytes());
//         kiley::bialignment::guided::global_guided(contig, &seq, &kiley_op, RADIUS, PARAMS).1
//     };
//     to_sam_ops(&kiley_op)
// }

// fn to_kiley_ops(cigar: &[sam::Op]) -> Vec<kiley::Op> {
//     cigar
//         .iter()
//         .flat_map(|op| match op {
//             sam::Op::Match(l) | sam::Op::Mismatch(l) | sam::Op::Align(l) => {
//                 std::iter::repeat(kiley::Op::Match).take(*l)
//             }
//             sam::Op::Deletion(l) => std::iter::repeat(kiley::Op::Del).take(*l),
//             sam::Op::Insertion(l) => std::iter::repeat(kiley::Op::Ins).take(*l),
//             _ => panic!("{:?}", op),
//         })
//         .collect()
// }

// fn to_sam_ops(k_ops: &[kiley::Op]) -> Vec<sam::Op> {
//     if k_ops.is_empty() {
//         return vec![];
//     }
//     fn is_the_same(op1: kiley::Op, op2: kiley::Op) -> bool {
//         match (op1, op2) {
//             (kiley::Op::Match, kiley::Op::Mismatch) | (kiley::Op::Mismatch, kiley::Op::Match) => {
//                 true
//             }
//             (x, y) if x == y => true,
//             _ => false,
//         }
//     }
//     assert!(!k_ops.is_empty());
//     let (mut current_op, mut len) = (k_ops[0], 1);
//     let mut ops = vec![];
//     for &op in k_ops.iter().skip(1) {
//         if is_the_same(op, current_op) {
//             len += 1;
//         } else {
//             match current_op {
//                 kiley::Op::Del => ops.push(sam::Op::Deletion(len)),
//                 kiley::Op::Ins => ops.push(sam::Op::Insertion(len)),
//                 kiley::Op::Mismatch | kiley::Op::Match => ops.push(sam::Op::Match(len)),
//             }
//             current_op = op;
//             len = 1;
//         }
//     }
//     match current_op {
//         kiley::Op::Del => ops.push(sam::Op::Deletion(len)),
//         kiley::Op::Ins => ops.push(sam::Op::Insertion(len)),
//         kiley::Op::Mismatch | kiley::Op::Match => ops.push(sam::Op::Match(len)),
//     }
//     ops
// }

fn to_op_runs(cigar: &[sam::Op]) -> Vec<(usize, char)> {
    cigar
        .iter()
        .filter_map(|op| match op {
            sam::Op::Match(l) | sam::Op::Mismatch(l) | sam::Op::Align(l) => Some((*l, 'M')),
            sam::Op::Deletion(l) => Some((*l, 'D')),
            _ => None,
        })
        .collect()
}

fn compress_small_matches(ops: &[(usize, char)], squish: usize) -> Vec<(usize, char)> {
    let ops = ops.iter().map(|&(len, op)| match op {
        'M' if len < squish => (len, 'D'),
        _ => (len, op),
    });
    compress(ops)
}

fn compress_small_deletions(ops: &[(usize, char)], squish: usize) -> Vec<(usize, char)> {
    let ops = ops.iter().map(|&(len, op)| match op {
        'D' if len < squish => (len, 'M'),
        _ => (len, op),
    });
    compress(ops)
}

fn compress<I: std::iter::Iterator<Item = (usize, char)>>(iter: I) -> Vec<(usize, char)> {
    let mut tidy_ops = vec![];
    let mut cum_len = 0;
    let mut prev = None;
    for (len, op) in iter {
        if prev.is_none() {
            prev = Some(op);
            cum_len = len;
        } else if Some(op) == prev {
            cum_len += len;
        } else {
            tidy_ops.push((cum_len, prev.unwrap()));
            cum_len = len;
            prev = Some(op);
        }
    }
    if let Some(op) = prev {
        tidy_ops.push((cum_len, op));
    }
    tidy_ops
}
use svg::node::element;
use svg::node::element::path::Data;

impl TidyAlignments {
    fn to_svg(&self, scale: &Scale, y_origin: usize) -> svg::node::element::Group {
        let y_position = y_origin + self.row * HEIGHT_PER_READS;
        let start_pos_bp = self.start;
        let end_pos_bp = self.end;
        // Add main line.
        let main_line = Data::new()
            .move_to((scale.map(start_pos_bp), y_position))
            .line_to((scale.map(end_pos_bp), y_position));
        let main_line = element::Path::new()
            .set("fill", "none")
            .set("stroke", "skyblue")
            .set("stroke-dasharray", 1)
            .set("stroke-width", DEL_STROKE)
            .set("stroke-opacity", DEL_OPACITY)
            .set("d", main_line);
        // Add clips.
        let clip_start = start_pos_bp.saturating_sub(self.clips.0);
        let head_clip = Data::new()
            .move_to((scale.map(clip_start), y_position))
            .line_to((scale.map(start_pos_bp), y_position));
        let head_clip = element::Path::new()
            .set("fill", "none")
            .set("stroke", "red")
            .set("stroke-width", DEL_STROKE)
            .set("stroke-opacity", CLIP_OPACITY)
            .set("d", head_clip);
        let clip_end = end_pos_bp + self.clips.1;
        let tail_clip = Data::new()
            .move_to((scale.map(end_pos_bp), y_position))
            .line_to((scale.map(clip_end), y_position));
        let tail_clip = element::Path::new()
            .set("fill", "none")
            .set("stroke", "red")
            .set("stroke-width", DEL_STROKE)
            .set("stroke-opacity", CLIP_OPACITY)
            .set("d", tail_clip);
        // Add matches.
        let mut match_rects = element::Group::new()
            .set("fill", "black")
            .set("stroke-width", 0);
        let mut rpos = start_pos_bp;
        for &(len, op) in self.operations.iter() {
            if op == 'M' {
                let width = scale.width(len);
                let rect = element::Rectangle::new()
                    .set("x", scale.map(rpos))
                    .set("width", width)
                    .set("y", y_position - MATCH_HEIGHT / 2)
                    .set("height", MATCH_HEIGHT);
                match_rects = match_rects.add(rect);
            }
            rpos += len;
        }
        let id = svg::node::Value::from(self.read_id.as_str());
        svg::node::element::Group::new()
            .set("id", id)
            .add(main_line)
            .add(tail_clip)
            .add(head_clip)
            .add(match_rects)
    }
}

const WIDTH: usize = 1500;
const HEIGHT_PER_READS: usize = 10;
const DEL_STROKE: usize = 3;
// Should be even.
const MATCH_HEIGHT: usize = 8;
const DEL_OPACITY: f64 = 0.7;
const CLIP_OPACITY: f64 = 0.3;
const SCALE_MARGIN: usize = 20;
const TOP_MARGIN: usize = 60;
const PILEUP_MARGIN: usize = SCALE_MARGIN + TOP_MARGIN;
const SIDE_MARGIN: usize = 50;
const LABEL_OFFSET: usize = 75;
const TICK_LEN: usize = 10;
const TICK: usize = 50_000;
const TICK_STROKE: usize = 2;
const SCALE_STROKE: usize = 5;
const SCALE_FONT_SIZE: usize = 25;
const TICK_FONT_SIZE: usize = 20;
fn flush_alignments(
    tidy_alignments: &[(String, Vec<TidyAlignments>)],
    selection: &Selection,
    args: &Args,
) {
    // for (id, alns) in tidy_alignments {
    //     eprintln!("==========={id}===============");
    //     for aln in alns.iter() {
    //         eprintln!("{aln:?}");
    //     }
    // }
    let scale = Scale::new(selection, args);
    let acc_heights: Vec<usize> = tidy_alignments
        .iter()
        .map(|(_, alns)| {
            let max = alns.iter().map(|x| x.row).max().unwrap_or(0);
            max * HEIGHT_PER_READS + PILEUP_MARGIN
        })
        .fold(vec![PILEUP_MARGIN], |mut acc, len| {
            acc.push(acc.last().unwrap() + len);
            acc
        });
    let height = *acc_heights.last().unwrap();
    let mut document = svg::Document::new()
        .set("viewBox", (0, 0, WIDTH, height))
        .set("width", WIDTH)
        .set("height", height);
    for (i, (id, alns)) in tidy_alignments.iter().enumerate() {
        let y_origin = acc_heights[i];
        let mut group = svg::node::element::Group::new();
        for aln in alns.iter() {
            group = group.add(aln.to_svg(&scale, y_origin));
        }
        group = group.add(scale.to_svg(id, y_origin));
        document = document.add(group);
    }
    svg::save(&args.output, &document).unwrap();
}

#[derive(Debug, Clone)]
struct Scale {
    // base pair.
    start: usize,
    end: usize,
    scale: f64,
    side_margin: f64,
}

impl Scale {
    fn to_svg(&self, id: &str, y_origin: usize) -> element::Group {
        let y_position = y_origin - SCALE_MARGIN;
        let id_node = svg::node::Value::from(id);
        let mut scale = svg::node::element::Group::new().set("id", id_node);
        let start = self.map(self.start);
        let end = self.map(self.end);
        let main_scale = Data::new()
            .move_to((start, y_position))
            .line_to((end, y_position));
        let main_scale = element::Path::new()
            .set("fill", "none")
            .set("stroke", "black")
            .set("stroke-width", SCALE_STROKE)
            .set("stroke-opacity", 1)
            .set("d", main_scale);
        scale = scale.add(main_scale);
        let label = svg::node::Text::new(id);
        let label = element::Text::new()
            .add(label)
            .set("x", start + LABEL_OFFSET as f64)
            .set("y", y_position - SCALE_STROKE - 5)
            .set("font-size", SCALE_FONT_SIZE);
        scale = scale.add(label);
        let tick_y_position = y_position - TICK_LEN;
        let tic_pos = (0..)
            .map(|i| i * TICK)
            .skip_while(|&pos| pos < self.start)
            .take_while(|&pos| pos <= self.end);
        for pos in tic_pos {
            let tick = Data::new()
                .move_to((self.map(pos), tick_y_position))
                .line_to((self.map(pos), y_position));
            let tick = element::Path::new()
                .set("fill", "none")
                .set("stroke", "black")
                .set("stroke-width", TICK_STROKE)
                .set("d", tick);
            scale = scale.add(tick);
            let label = format!("{} Kbp", pos / 1_000);
            let label = svg::node::Text::new(label);
            let label = element::Text::new()
                .add(label)
                .set("x", self.map(pos))
                .set("y", y_position - SCALE_STROKE - 5)
                .set("font-size", TICK_FONT_SIZE);
            scale = scale.add(label);
        }
        scale
    }
    fn width(&self, width_in_bp: usize) -> f64 {
        tidy(self.scale * width_in_bp as f64)
    }
    fn map(&self, pos_in_bp: usize) -> f64 {
        let pos_in_bp = pos_in_bp.max(self.start).min(self.end);
        let pos_in_pix = self.side_margin as f64 + (pos_in_bp - self.start) as f64 * self.scale;
        tidy(pos_in_pix)
    }
    fn new(selection: &Selection, _args: &Args) -> Self {
        let start = selection.selects.iter().map(|x| x.1).min().unwrap();
        let end = selection.selects.iter().map(|x| x.2).max().unwrap();
        let length = end - start;
        let width = (WIDTH - SIDE_MARGIN * 2) as f64;
        let scale = width / length as f64;
        let side_margin = SIDE_MARGIN as f64;
        Self {
            start,
            end,
            scale,
            side_margin,
        }
    }
}

fn tidy(x: f64) -> f64 {
    (x * 1000f64).ceil() / 1000f64
}

use bio_utils::sam::Sam;
use std::collections::HashMap;
type Samfile = (HashMap<String, usize>, Vec<Sam>);
fn parse_sam_file<I: std::iter::Iterator<Item = String>>(lines: I) -> Samfile {
    let mut contig_length = HashMap::new();
    let mut alignments = vec![];
    for line in lines {
        if line.starts_with("@SQ") {
            let line: Vec<_> = line.split('\t').collect();
            assert!(line[1].starts_with("SN:"));
            let (_, tig_name) = line[1].split_once(':').unwrap();
            assert!(line[2].starts_with("LN:"));
            let (_, length) = line[2].split_once(':').unwrap();
            let length: usize = length.parse().unwrap();
            contig_length.insert(tig_name.to_string(), length);
        } else if !line.starts_with('@') {
            let sam = Sam::new(&line).unwrap();
            alignments.push(sam);
        }
    }
    (contig_length, alignments)
}
