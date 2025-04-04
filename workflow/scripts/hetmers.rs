use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use sha2::{Digest, Sha256};
use std::env;
use std::process;

fn get_hetmers(input: &str, alleles: &[usize], minimum: usize, output_prefix: &str) {
    println!("Loading k-mer count file {}...", input);
    let file = File::open(input).expect("Unable to open file");
    let reader = BufReader::new(file);

    let mut seqs = Vec::new();
    let mut counts = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            continue;
        }
        let seq = parts[0].to_string();
        let count: usize = parts[1].parse().expect("Invalid count value");
        if count >= minimum {
            seqs.push(seq);
            counts.push(count);
        }
    }

    let k = seqs[0].len();
    println!("k is {}", k);
    let k_half = k / 2;

    let no_center_ks: Vec<String> = seqs.iter()
        .map(|s| format!("{}{}", &s[..k_half], &s[k_half + 1..]))
        .collect();

    println!("Reverse complementing...");
    let complement = |c: char| match c {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c
    };
    let no_center_ks_rev: Vec<String> = no_center_ks.iter()
        .map(|s| s.chars().rev().map(complement).collect())
        .collect();

    println!("Hashing...");
    let hash_fn = |s: &String| {
        let mut hasher = Sha256::new();
        hasher.update(s.as_bytes());
        u64::from_be_bytes(hasher.finalize()[..8].try_into().unwrap())
    };
    
    let no_center_ks_hash: Vec<u64> = no_center_ks.iter().map(hash_fn).collect();
    let no_center_ks_rev_hash: Vec<u64> = no_center_ks_rev.iter().map(hash_fn).collect();

    println!("Getting the minimum hash...");
    let min_hashes: Vec<u64> = no_center_ks_hash.iter()
        .zip(no_center_ks_rev_hash.iter())
        .map(|(x, y)| *x.min(y))
        .collect();
    
    println!("Grouping unique hashes into a dictionary...");
    let mut d: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, num) in min_hashes.iter().enumerate() {
        d.entry(*num).or_insert_with(Vec::new).push(i);
    }

    println!("Filtering hashes...");
    let ans: HashMap<u64, Vec<usize>> = d.into_iter()
        .filter(|(_, v)| alleles.contains(&v.len()))
        .collect();

    println!("Extracting counts and sequences...");
    let hetmers: Vec<String> = ans.values()
        .map(|indices| indices.iter().map(|&i| seqs[i].clone()).collect::<Vec<String>>().join(","))
        .collect();
    let hetmer_counts: Vec<String> = ans.values()
        .map(|indices| indices.iter().map(|&i| counts[i].to_string()).collect::<Vec<String>>().join(","))
        .collect();

    println!("Saving results...");
    let mut seq_file = File::create(format!("{}_seqs.csv", output_prefix)).expect("Unable to create file");
    writeln!(seq_file, "{}", hetmers.join("\n")).expect("Unable to write to file");
    
    let mut count_file = File::create(format!("{}_counts.csv", output_prefix)).expect("Unable to create file");
    writeln!(count_file, "{}", hetmer_counts.join("\n")).expect("Unable to write to file");
    
    println!("Done! :D");
}


fn main() {
    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();

    // Check if the correct number of arguments is provided
    if args.len() < 5 {
        eprintln!("Usage: {} <input> <alleles> <minimum> <output_prefix>", args[0]);
        process::exit(1);
    }

    let input = &args[1];  // Input file path
    let alleles: Vec<usize> = args[2]
        .split(',')
        .filter_map(|s| s.parse().ok())
        .collect();  // Parse allele numbers
    let minimum: usize = args[3]
        .parse()
        .expect("Minimum should be a positive integer");
    let output_prefix = &args[4];  // Output file prefix

    // Call the function
    get_hetmers(input, &alleles, minimum, output_prefix);
}
