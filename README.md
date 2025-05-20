# poolseq-kmers

Author: Miles Roberts

Simulation workflow to investigate the utility of k-mers, het-mers, and k-unitigs for pool-seq data analysis built with snakemake (v 7.25.0)

## Table of Contents

## Overview 

## Inputs

See config/README.md for a complete description of workflow inputs. In short, you need two files:

* config/config.yaml:

* config/parameters.tsv:

## Outputs

For each simulation, this workflow outputs:

* SNP calls from Varscan

* SNP calls from PoolSNP

* SNP calls from discosnp

* Het-mers from smudgeplot

* Het-mers from hetmers

## Usage

Examples commands are in resources/01_snakemake.sh

### Run whole workflow with conda envs on slurm cluster

### Run workflow in batches

### Run workflow with singularity

Need to pass `--use-singularity` to snakemake and also your snakemake working directory with `--singularity-args "--bind <SNAKEMAKE_WORKING_DIRECTORY>"`

```
snakemake --cores 1 --use-singularity --singularity-args "--bind ~/Josephs_Lab_Projects/poolseq-kmers/workflow"
```

## Notes

### Building docker container

https://github.com/snakemake/snakemake/issues/2602

Create docker file with `snakemake --containerize > Dockerfile`. Copy Dockerfile and envs to same directory. Then run these commands from a computer with docker:

```
sudo docker build -t poolseq-kmers .
sudo docker login -u milesroberts
sudo docker tag poolseq-kmers milesroberts/poolseq-kmers
sudo docker push milesroberts/poolseq-kmers
```

## To do

### higher priority

- [x] add option to vary sequencing machine

- [x] add rule to remove regions from references

- [x] add discosnp for comparison

- [x] add poolsnp: https://github.com/capoony/PoolSNP

- [x] add software to call unitigs then align them back to reference genome 

- [x] compare snp and genome-wide diversity estimates to ground truth

- [x] compare snp and genome-wide fst estimates to ground truth

- [x] simulate reference bias by masking portions of reference genome (how much to mask, which individuals to mask)

- [x] add hetmers binary to workflow

- [x] generalize bayes theorem to work with minimum kmer counts > 1

- [x] generate genomes that follow power-law distributions for k-mer counts

- [x] put rule outputs in separate folders

- [x] add script to model QTLs in a mapping population

- [x] mark kmc output files as temp files

- [x] add more parameters to config.yaml

- [x] parallelize discosnp for one population

- [x] parallelize discosnp for two populations

- [x] add bulk segregant analysis simulation for fst

- [x] output hetmers to their own directory

- [x] figure out way to run workflow in batches

- [x] remove need for ref.txt in config/

- [x] write workflow schema for config.yaml

- [x] write workflow schema for parameters.tsv

- [x] resolve workflow lints `snakemake --lint`

- [x] figure out singularity

- [ ] upgrade to latest snakemake version

- [ ] integration tests

- [ ] github actions

- [ ] write hetmers binary to calculate fst

- [ ] generalize bayes theorem to negative binomial distribution

### lower priority

- [ ] add in ploidyfrost

- [ ] try adding kmer2snp

- [ ] add purifying selection simulation - does this also create unitigs?

- [ ] update to smudgeplot >0.3.0, once we're able to get k-mer sequences again

- [ ] figure out how to use unpaired reads in varscan 
