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


