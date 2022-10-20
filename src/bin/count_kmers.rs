use std::path::PathBuf;

use clap::Parser;

/// Dump k-mer histograms in TSV format. It is the canonicalized counts, i.e., forward strands and reverse strands are merged together.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Genomes<FASTA>.
    #[arg(short, long)]
    genomes: PathBuf,
    /// Minimum k-mer.
    #[arg(short, long, default_value_t = 10)]
    min_k_mer: usize,
    /// Maximum K-mer.
    #[arg(short, long, default_value_t = 20)]
    max_k_mer: usize,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let genomes: Vec<_> = bio_utils::fasta::parse_into_vec(&args.genomes)?;
    println!("K\tOccInGenome\tNumOfKmer");
    for k in args.min_k_mer..=args.max_k_mer {
        let histogram = counts(&genomes, k);
        for (count, freq) in histogram {
            println!("{k}\t{count}\t{freq}");
        }
    }
    Ok(())
}

use std::collections::HashMap;
fn counts(genomes: &[bio_utils::fasta::Record], k: usize) -> Vec<(u32, u32)> {
    let counts = mito_check::count_kmers(genomes, k);
    let mut freq_in_genomes: HashMap<_, u32> = HashMap::new();
    for &count in counts.values() {
        *freq_in_genomes.entry(count).or_default() += 1;
    }
    let mut freq_in_genomes: Vec<_> = freq_in_genomes.into_iter().collect();
    freq_in_genomes.sort_unstable_by_key(|x| x.0);
    freq_in_genomes
}
