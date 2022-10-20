use std::path::PathBuf;

use clap::Parser;

/// Enumerate repetitive kmers in the genomes in TSV.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Genomes<FASTA>.
    #[arg(short, long)]
    genomes: PathBuf,
    /// K for the k-mer.
    #[arg(short, long, default_value_t = 13)]
    kmer: usize,
    /// Threshold for determine repetitiveness.
    #[arg(short, long, default_value_t = 15)]
    threshold: u32,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let genomes = bio_utils::fasta::parse_into_vec(&args.genomes)?;
    let repetitive_kmers = {
        let mut counts = mito_check::count_kmers(&genomes, args.kmer);
        counts.retain(|_, &mut count| args.threshold <= count);
        counts
    };
    println!("ID\tStart\tEnd\tCount\tSeq");
    for genome in genomes.iter() {
        let mut match_region = None;
        for (i, kmer) in genome.seq().windows(args.kmer).enumerate() {
            let idx = mito_check::to_idx(kmer);
            let count = match repetitive_kmers.get(&idx) {
                Some(count) => *count,
                None => continue,
            };
            match match_region {
                None => match_region = Some((i, i, count)),
                Some((start, end, total)) if end + 1 == i => {
                    match_region = Some((start, i, total + count));
                }
                Some((start, end, total)) => {
                    let seq = std::str::from_utf8(&genome.seq()[start..end + args.kmer]).unwrap();
                    let ave_count = total / (end - start + 1) as u32;
                    println!("{}\t{start}\t{end}\t{ave_count}\t{seq}", genome.id(),);
                    match_region = Some((i, i, count));
                }
            }
        }
    }
    Ok(())
}
