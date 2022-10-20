use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The pileup.
    #[arg(short, long)]
    pileup: PathBuf,
    /// Minimum support for variants.
    #[arg(short, long, default_value_t = 2)]
    min_support: usize,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    let pileup = bio_utils::fasta::parse_into_vec(&args.pileup).unwrap();
    let num_var = num_supported_variants(&pileup, args.min_support);
    if 0 < num_var {
        let sid = pileup[0].id();
        println!("{sid}\t{num_var}");
        std::process::ExitCode::SUCCESS
    } else {
        std::process::ExitCode::FAILURE
    }
}

use std::collections::HashMap;
fn num_supported_variants(pileup: &[bio_utils::fasta::Record], min: usize) -> usize {
    let len = pileup[0].seq().len();
    let mut summary: Vec<HashMap<_, usize>> = vec![HashMap::new(); len];
    for record in pileup.iter() {
        for (i, &base) in record.seq().iter().enumerate() {
            if base != b'-' {
                *summary[i].entry(base).or_default() += 1;
            }
        }
    }
    summary
        .iter()
        .filter_map(|hm| {
            let mut values: Vec<_> = hm.values().copied().collect();
            values.sort();
            values.pop();
            values.pop()
        })
        .filter(|&count| min <= count)
        .count()
}
