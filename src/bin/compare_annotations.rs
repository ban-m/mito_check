use clap::Parser;
use std::path::PathBuf;

/// Compare break points and repetitive regions statistically.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Genomes
    #[arg(short, long)]
    genomes: PathBuf,
    /// Repeat Annotation.
    #[arg(short, long)]
    repeat_annotation: PathBuf,
    /// Break point annotation.
    #[arg(short, long)]
    break_point_annotation: PathBuf,
    /// Output prefix.
    #[arg(short, long)]
    output_prefix: PathBuf,
    /// Seed
    #[arg(short, long, default_value_t = 348290)]
    seed: u64,
}

use std::collections::HashMap;
const SAMPLE_NUM: u64 = 100_000;
const MINI_SAMPLE: u64 = 100;
use std::io::*;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let genomes = bio_utils::fasta::parse_into_vec(&args.genomes)?;
    let mut repeat_annotation: HashMap<_, Vec<_>> = HashMap::new();
    for line in std::fs::File::open(&args.repeat_annotation)
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .skip(1)
    // Skip header.
    {
        let line: Vec<_> = line.split('\t').collect();
        let start: usize = line[1].parse().unwrap();
        let end: usize = line[2].parse().unwrap();
        assert!(start <= end);
        repeat_annotation
            .entry(line[0].to_string())
            .or_default()
            .push((start, end));
    }
    repeat_annotation
        .values_mut()
        .for_each(|xs| xs.sort_by_key(|x| x.0));
    let break_points: Vec<_> = std::fs::File::open(&args.break_point_annotation)
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .skip(1) // Skip header.
        .map(|line| {
            let line: Vec<_> = line.split('\t').collect();
            let start: usize = line[1].parse().unwrap();
            let end: usize = line[2].parse().unwrap();
            (line[0].to_string(), start, end)
        })
        .collect();
    let nearest_repeat = get_nearest_repeats(&repeat_annotation, &break_points);
    {
        let mut outfile = args.output_prefix.clone();
        outfile.push("dist_to_nearest_repeats.tsv");
        let mut wtr = std::fs::File::create(outfile).map(BufWriter::new)?;
        writeln!(&mut wtr, "ID\tDistance\tIsNull")?;
        for (contig, distance) in nearest_repeat.iter() {
            writeln!(&mut wtr, "{contig}\t{distance}\tObs")?;
        }
        let null_distr: Vec<_> = (0..MINI_SAMPLE)
            .flat_map(|i| {
                let seed = args.seed + i + SAMPLE_NUM;
                sample_nearest_repeats(&repeat_annotation, &break_points, &genomes, seed)
            })
            .collect();
        for (contig, distance) in null_distr {
            writeln!(&mut wtr, "{contig}\t{distance}\tSampled")?;
        }
    }
    let average_distance = {
        let total: usize = nearest_repeat.iter().map(|x| x.1).sum();
        total / nearest_repeat.len()
    };
    let null_distr: Vec<_> = (0..SAMPLE_NUM)
        .map(|i| {
            let seed = args.seed + i as u64;
            let null_distr =
                sample_nearest_repeats(&repeat_annotation, &break_points, &genomes, seed);
            let total: usize = null_distr.iter().map(|x| x.1).sum();
            total / null_distr.len()
        })
        .collect();
    {
        let mut outfile = args.output_prefix;
        outfile.push("average_dist_to_repeats.tsv");
        let mut wtr = std::fs::File::create(outfile).map(BufWriter::new)?;
        writeln!(&mut wtr, "Distance\tIsNull")?;
        writeln!(&mut wtr, "{average_distance}\tObs")?;
        for distance in null_distr.iter() {
            writeln!(&mut wtr, "{distance}\tSampled")?;
        }
    }
    let above_average = null_distr
        .iter()
        .filter(|&&ave| ave < average_distance)
        .count();
    println!("{}", (above_average as f64) / (SAMPLE_NUM as f64));
    Ok(())
}

type RepeatAnnots = HashMap<String, Vec<(usize, usize)>>;
type BreakAnnot = (String, usize, usize);
fn get_nearest_repeats(
    repeat_annotations: &RepeatAnnots,
    break_points: &[BreakAnnot],
) -> Vec<(String, usize)> {
    let mut dist_to_nearest_repeats = vec![];
    for (contig, start, end) in break_points.iter() {
        let repeats = repeat_annotations.get(contig).unwrap();
        dist_to_nearest_repeats.push((contig.to_string(), find_nearest(repeats, *start)));
        dist_to_nearest_repeats.push((contig.to_string(), find_nearest(repeats, *end)));
    }
    dist_to_nearest_repeats
}

fn find_nearest(repeats: &[(usize, usize)], position: usize) -> usize {
    repeats
        .iter()
        .map(|&(start, end)| {
            if start < position && position < end {
                0
            } else {
                abs(start, position).min(abs(end, position))
            }
        })
        .min()
        .unwrap()
}

fn abs(x: usize, y: usize) -> usize {
    x.max(y) - x.min(y)
}

use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
fn sample_nearest_repeats(
    repeat_annotations: &RepeatAnnots,
    break_points: &[BreakAnnot],
    genomes: &[bio_utils::fasta::Record],
    seed: u64,
) -> Vec<(String, usize)> {
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let mut null_distr = Vec::with_capacity(break_points.len() * 2);
    for (contig, _, _) in break_points.iter() {
        let annotation = repeat_annotations.get(contig).unwrap();
        let genome = genomes.iter().find(|g| g.id() == contig).unwrap();
        let len = genome.seq().len();
        let dist = find_nearest(annotation, rng.gen_range(0..len));
        null_distr.push((contig.to_string(), dist));
        let dist = find_nearest(annotation, rng.gen_range(0..len));
        null_distr.push((contig.to_string(), dist));
    }
    null_distr
}
