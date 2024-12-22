# poolseq-kmers

trying out approaches to analyzing k-mers in poolseq data

## IDEAS

I'm thinking that I can identify pairs of k-mers that differ at their central bp as putative snps, then use the relative coverage between the k-mers to estimate minor allele frequencies.

If I can get a site frequency spectrum (could compare k-mers to ancestral/outgroup k-mers to determine which is derived and which is ancestral)

What about when individuals do not contribute equally to the genome pool?

What about when the genome sequence is not repetitive? I could vary the shannon entropy of the sequences I use

## workflow

* generate ancestral sequence

* neutral forward-time simulation in slim

* generate sequencing reads with insilicoseq

* count k-mers with kmc

* use smudgeplot to get k-mers that differ at central base pair

* calculate minor allele frequencies based on minor k-mer coverage

* map k-mer coverages back to actual sequences using seqkit

## parameters

* type of ancestral sequence

* forward simulation parameters

* number of individuals per pool, error model, coverage, variation in individual contribution

* coverage cutoff for k-mers

## notes

Dang, smudgeplot is going through a big update right now, which is just my luck. 

### 2024-12-21

Trying to debug the workflow, but it's very complex with lots of paths. I'll just debug chunks of the workflow at a time

A snippet for debugging:

```
snakemake --cluster "sbatch --time={resources.time} --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb_per_cpu} --partition=josephsnodes --account=josephsnodes" --resources load=1 --jobs 950 --cores 950 --use-conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --keep-incomplete calls_2600_p1.tsv kmerpairs_2600_p1_coverages.tsv kmerpairs_2600_p1_sequences.tsv discoRes_ad_2600_p1.txt slim_allele_freqs_2600_p1.txt 2600_p1_poolsnp_output.vcf.gz calls_2600_p2.tsv kmerpairs_2600_p2_coverages.tsv kmerpairs_2600_p2_sequences.tsv discoRes_ad_2600_p2.txt slim_allele_freqs_2600_p2.txt 2600_p2_poolsnp_output.vcf.gz
```

## to do

- [x] add option to vary sequencing machine

- [x] add rule to remove regions from references

- [x] add discosnp for comparison

- [x] add poolsnp: https://github.com/capoony/PoolSNP

- [ ] update to smudgeplot >0.3.0, once we're able to get k-mer sequences again

- [ ] figure out how to use unpaired reads in varscan 