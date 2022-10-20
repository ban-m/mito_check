use bio_utils::paf::PAF;
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Alignments<PAF with cigar tags>.
    #[arg(short, long)]
    alignments: PathBuf,
    /// Prefix of the output path.
    #[arg(short, long)]
    prefix: PathBuf,
    /// Upper limit of the insertion/deletion
    #[arg(short, long, default_value_t = 1000)]
    upper_size_of_indel: usize,
}
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    use std::io::*;
    let mut diff = vec![];
    let records: Vec<_> = std::fs::File::open(&args.alignments)
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
        .filter_map(|l| bio_utils::paf::PAF::new(&l))
        .collect();
    for record in records.iter() {
        let mut qpos = 0;
        let mut rpos = record.tstart;
        let qname = &record.qname;
        let rname = &record.tname;
        let cigar = bio_utils::sam::parse_cigar_string(record.get_tag("cg").unwrap().1);
        for op in cigar {
            use bio_utils::sam::Op;
            match op {
                Op::Match(l) | Op::Align(l) => {
                    qpos += l;
                    rpos += l;
                }
                Op::Insertion(l) => {
                    if record.relstrand {
                        let qpos = record.qstart + qpos;
                        diff.push((qname, qpos, rname, rpos, "Ins", l));
                    } else {
                        let qpos = record.qend - qpos;
                        diff.push((qname, qpos, rname, rpos, "Ins", l));
                    }
                    qpos += l;
                }
                Op::Deletion(l) => {
                    if record.relstrand {
                        let qpos = record.qstart + qpos;
                        diff.push((qname, qpos, rname, rpos, "Del", l));
                    } else {
                        let qpos = record.qend - qpos;
                        diff.push((qname, qpos, rname, rpos, "Del", l));
                    }
                    rpos += l;
                }
                Op::Mismatch(l) => {
                    if record.relstrand {
                        let qpos = record.qstart + qpos;
                        diff.push((qname, qpos, rname, rpos, "Mism", l));
                    } else {
                        let qpos = record.qend - qpos;
                        diff.push((qname, qpos, rname, rpos, "Mism", l));
                    }
                    rpos += l;
                    qpos += l;
                }
                _ => {}
            }
        }
    }
    {
        let mut filename = args.prefix.clone();
        filename.push("diff.tsv");
        let mut wtr = std::fs::File::create(&filename).map(BufWriter::new)?;
        writeln!(&mut wtr, "Query\tQpos\tRefr\tRpos\tType\tSize")?;
        for (qname, qpos, rname, rpos, t, size) in diff.iter() {
            writeln!(&mut wtr, "{qname}\t{qpos}\t{rname}\t{rpos}\t{t}\t{size}")?;
        }
    }
    let error_rates = get_error_profile(&records, &diff, args.upper_size_of_indel);
    {
        let mut filename = args.prefix;
        filename.push("diff_rate.tsv");
        let mut wtr = std::fs::File::create(&filename).map(BufWriter::new)?;
        for (cid, ins, del, mism, len) in error_rates.iter() {
            writeln!(&mut wtr, "{cid}\t{ins}\t{del}\t{mism}\t{len}")?;
        }
    }
    let no_error: f64 = error_rates
        .iter()
        .map(|(_, _, _, mism, _)| (1f64 - mism))
        .product();
    let single_error: f64 = error_rates
        .iter()
        .enumerate()
        .map(|(i, &(_, _, _, mism, _))| {
            let no_error: f64 = error_rates
                .iter()
                .enumerate()
                .filter_map(|(j, (_, _, _, m, _))| (i != j).then_some(1f64 - *m))
                .product();
            no_error * mism
        })
        .sum();
    let two_or_more_error = 1f64 - no_error - single_error;
    println!("{two_or_more_error}");
    Ok(())
}

type ErrorProfile = (String, f64, f64, f64, usize);
fn get_error_profile(
    alns: &[PAF],
    diff: &[(&String, usize, &String, usize, &str, usize)],
    max_indel: usize,
) -> Vec<ErrorProfile> {
    use std::collections::HashMap;
    let strain_to_length: HashMap<_, _> = {
        let contig_to_length: HashMap<_, _> = alns
            .iter()
            .map(|paf| (paf.qname.to_string(), paf.qlen))
            .collect();
        let mut strain_to_length: HashMap<_, usize> = HashMap::new();
        for (contig, length) in contig_to_length {
            let (strain, _) = contig.split_once('_').unwrap();
            *strain_to_length.entry(strain.to_string()).or_default() += length;
        }
        strain_to_length
    };
    let mut strain_diffs: HashMap<_, _> = strain_to_length
        .keys()
        .map(|st| (st.clone(), (0, 0, 0)))
        .collect();
    for (qname, _, _, _, err, size) in diff.iter().filter(|x| x.5 < max_indel) {
        let (strain, _) = qname.split_once('_').unwrap();
        let slot = strain_diffs.get_mut(strain).unwrap();
        match *err {
            "Ins" => slot.0 += size,
            "Del" => slot.1 += size,
            "Mism" => slot.2 += size,
            _ => panic!(),
        }
    }
    strain_diffs
        .into_iter()
        .map(|(strain, (ins, del, mism))| {
            let length = strain_to_length[&strain];
            let ins = ins as f64 / length as f64;
            let del = del as f64 / length as f64;
            let mism = mism as f64 / length as f64;
            (strain, ins, del, mism, length)
        })
        .collect()
}
