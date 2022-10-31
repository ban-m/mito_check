use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;

/// Summarize assemblies and annotations.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Assembled contigs in FASTA format.
    #[arg(short, long)]
    contigs: PathBuf,
    /// Annotation in GFF3 format (*Reference*).
    #[arg(short, long)]
    gff: PathBuf,
    /// Alignments in last MAF format.
    #[arg(short, long)]
    maf: PathBuf,
    /// Output variant information into.
    #[arg(short, long)]
    output: PathBuf,
    /// Filter out alignments below this size.
    #[arg(short, long, default_value_t = 2000)]
    min_aln_size: u64,
}

use std::io::prelude::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let contigs = bio_utils::fasta::parse_into_vec(&args.contigs)?;
    let maf = {
        let mut maf = bio_utils::maf::parse_into_vec(&args.maf)?;
        maf.retain(|record| {
            let sequences = record.sequence();
            let refr = &sequences[0];
            let query = &sequences[1];
            args.min_aln_size < refr.length().min(query.length())
        });
        maf
    };
    let genome_size: usize = contigs.iter().map(|x| x.seq().len()).sum();
    let num_contigs = contigs.len();
    let mut variants: Vec<_> = vec![];
    for record in maf.iter() {
        let sequences = record.sequence();
        let refr = &sequences[0];
        let refname = refr.name();
        let strand = refr.is_forward();
        let length = refr.length();
        let start = refr.start();
        let query = &sequences[1];
        assert_eq!(refr.text().len(), query.text().len());
        let mut rpos = 0;
        for (&r, &q) in std::iter::zip(refr.text(), query.text()) {
            let (r, q) = (r.to_ascii_uppercase(), q.to_ascii_uppercase());
            if r == b'-' || q == b'-' || r != q {
                let var_type = if r == b'-' {
                    "Ins"
                } else if q == b'-' {
                    "Del"
                } else {
                    "Subs"
                };
                let position = match strand {
                    true => start + rpos,
                    false => length - start - rpos,
                };
                variants.push((refname, position, var_type));
            }
            rpos += (r != b'-') as u64;
        }
    }
    let variations = variants.len();
    let (cov_refr, cov_contigs) = get_coverage(&maf);
    println!("{genome_size}\t{num_contigs}\t{variations}\t{cov_refr}\t{cov_contigs}");
    // GFF
    let gff: Vec<_> = std::fs::File::open(&args.gff)
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
        .collect();
    let exons: Vec<Vec<_>> = gff
        .iter()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let line: Vec<_> = line.split('\t').collect();
            (line[2] == "exon").then_some(line)
        })
        .collect();
    let mut wtr = std::fs::File::create(args.output).map(BufWriter::new)?;
    for (rname, position, var_type) in variants {
        let hit_exon: Vec<_> = exons
            .iter()
            .filter(|exon| {
                let start: u64 = exon[3].parse().unwrap();
                let end: u64 = exon[4].parse().unwrap();
                exon[0] == rname && start < position && position <= end
            })
            .collect();
        for exon in hit_exon.iter() {
            writeln!(
                &mut wtr,
                "{rname}\t{position}\t{var_type}\t{}",
                exon.join("\t")
            )?;
        }
        if hit_exon.is_empty() {
            writeln!(&mut wtr, "{rname}\t{position}\t{var_type}\tNonExon")?;
        }
    }
    Ok(())
}

use std::collections::HashMap;

fn get_coverage(mafs: &[bio_utils::maf::Record]) -> (f64, f64) {
    let (refr_sizes, contig_sizes): (HashMap<_, _>, HashMap<_, _>) = mafs
        .iter()
        .map(|record| {
            let sequences = record.sequence();
            let refr = &sequences[0];
            let contig = &sequences[1];
            let refr = (refr.name().to_string(), refr.src_size());
            let contig = (contig.name().to_string(), contig.src_size());
            (refr, contig)
        })
        .unzip();
    let mut refr_covered: HashMap<_, _> = refr_sizes.keys().cloned().map(|x| (x, 0)).collect();
    let mut contig_covered: HashMap<_, _> = contig_sizes.keys().cloned().map(|x| (x, 0)).collect();
    for record in mafs {
        let sequences = record.sequence();
        let refr = &sequences[0];
        let contig = &sequences[1];
        *refr_covered.get_mut(refr.name()).unwrap() += refr.length();
        *contig_covered.get_mut(contig.name()).unwrap() += contig.length();
    }
    let refr_cov: u64 = refr_covered.values().sum();
    let refr_len: u64 = refr_sizes.values().sum();
    let contig_cov: u64 = contig_covered.values().sum();
    let contig_len: u64 = contig_sizes.values().sum();
    (
        refr_cov as f64 / refr_len as f64,
        contig_cov as f64 / contig_len as f64,
    )
}
