use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let stdin = std::io::stdin();
    let maf = bio_utils::maf::Reader::new(BufReader::new(stdin.lock()));
    let arg: Vec<_> = std::env::args().collect();
    let len: u64 = arg[1].parse().unwrap();
    let (mut ins, mut del, mut subs, mut aln) = (0, 0, 0, 0);
    for record in maf.records().filter_map(|x| x.ok()) {
        let sequences = record.sequence();
        let refr = &sequences[0];
        let query = &sequences[1];
        if refr.length() < len || query.length() < len {
            continue;
        }
        assert_eq!(refr.text().len(), query.text().len());
        for (&r, &q) in std::iter::zip(refr.text(), query.text()) {
            ins += (r == b'-') as usize;
            del += (q == b'-') as usize;
            subs += (r.to_ascii_uppercase() != q.to_ascii_uppercase()) as usize;
        }
        aln += ((refr.length() as f64).sqrt() * (query.length() as f64).sqrt()).ceil() as usize;
    }
    //   println!("INS\tDEL\tSUBS\tTOTAL\tALN");
    println!("{ins}\t{del}\t{subs}\t{}\t{aln}", ins + del + subs);
    Ok(())
}
