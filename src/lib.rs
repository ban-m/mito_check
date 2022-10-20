use std::collections::HashMap;
pub fn count_kmers(genomes: &[bio_utils::fasta::Record], k: usize) -> HashMap<u64, u32> {
    let mut counts: HashMap<_, u32> = HashMap::new();
    for genome in genomes.iter() {
        for kmer in genome.seq().windows(k) {
            let kmer = to_idx(kmer);
            *counts.entry(kmer).or_default() += 1;
        }
    }
    counts
}

pub fn back_to_seq(kmer: u64, k: usize) -> Vec<u8> {
    (0..k)
        .map(|offset| {
            let idx = (kmer >> (2 * offset)) & 0b11;
            b"ACGT"[idx as usize]
        })
        .rev()
        .collect()
}

pub fn to_idx(w: &[u8]) -> u64 {
    // Determine if this k-mer is canonical.
    let is_canonical = {
        let mut idx = 0;
        while idx < w.len() / 2
            && w[idx].to_ascii_uppercase() == w[w.len() - idx - 1].to_ascii_uppercase()
        {
            idx += 1;
        }
        w[idx] <= w[w.len() - idx - 1]
    };
    if is_canonical {
        w.iter()
            .fold(0, |cum, &x| (cum << 2) | BASE2BIT[x as usize])
    } else {
        w.iter()
            .rev()
            .fold(0, |cum, &x| (cum << 2) | BASE2BITCMP[x as usize])
    }
}

const BASE2BITCMP: [u64; 256] = base2bitcmp();
const BASE2BIT: [u64; 256] = base2bit();

const fn base2bitcmp() -> [u64; 256] {
    let mut slots = [0; 256];
    slots[b'A' as usize] = 3;
    slots[b'a' as usize] = 3;
    slots[b'C' as usize] = 2;
    slots[b'c' as usize] = 2;
    slots[b'G' as usize] = 1;
    slots[b'g' as usize] = 1;
    slots
}

const fn base2bit() -> [u64; 256] {
    let mut slots = [0; 256];
    slots[b'C' as usize] = 1;
    slots[b'c' as usize] = 1;
    slots[b'G' as usize] = 2;
    slots[b'g' as usize] = 2;
    slots[b'T' as usize] = 3;
    slots[b't' as usize] = 3;
    slots
}
