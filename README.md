# poolseq-kmers

[![Super-Linter](https://github.com/milesroberts-123/poolseq-kmers/actions/workflows/linter.yml/badge.svg)](https://github.com/marketplace/actions/super-linter) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io) [![GitHub actions status](https://github.com/milesroberts-123/poolseq-kmers/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests) [![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/) [![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>)

Simulation workflow to investigate the utility of k-mers, het-mers, and k-unitigs for pool-seq data analysis built with snakemake (v 9.3.3) and run with the snakemake slurm plugin (v 1.3.6).

Author: Miles Roberts

## Table of Contents

[Setup](#setup)

[Inputs](#inputs)

[Outputs](#outputs)

[Usage](#usage)

## Setup

1. Install conda/mamba

2. Download the workflow from github

3. Create a mamba environment with snakemake and any plugins you need

This is how to make new mamba environment named snakemake with snakemake and the slurm plugin installed. If you are not running the workflow on a SLURM cluster, you can install a different pluggin

```sh
conda create -y -f snakemake-mamba-env.yaml

conda activate snakemake
```

4. Check the snakemake profile for the proper executer. The default profile runs snakemake on a slurm cluster (`workflow/profiles/default/config.yaml`), but you should still change the slurm account, slurm partition, and default resources to match your system.

## Inputs

See config/README.md for a complete description of workflow inputs. In short, There are three main branches to the workflow `all_histos` which will generate k-mer histograms, `all_sims` which will run SLiM simulations, and `all_empirical` which will test the performance of freqk on provided emprical datasets. You can run all three branches at once with the target rule `all` (default). For any run, you need `config/config.yaml` and `config/parameters.tsv`.

## Outputs

The workflow outputs vary by branch

### all (default)

Runs all three branches of the workflow.

### all histos

### all sims

### all empirical

## Usage

There are three main branches to the workflow `all_histos` which will generate k-mer histograms, `all_sims` which will run SLiM simulations, and `all_empirical` which will test the performance of freqk on provided emprical datasets. You can run all three branches at once with the target rule `all` (default). 

In addition, the workflow can be run with conda environments or with docker (recommended).

If you're doing lots of simulations, then run snakemake in batches (example below)

Example commands are given below.

### Run whole workflow with conda envs on slurm cluster

`snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going`

### Run simulations in batches with conda envs on slurm cluster

```sh
batches=50
for i in {1..$batches}
do
  snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --keep-going --batch all_sims=$i/$batches all_sims
done
```

### Run workflow with singularity on a slurm cluster

Instead of downloading and building all of the conda environments, you can just download a container with all of the conda environments pre-installed.

Need to pass `--use-singularity` to snakemake and also your snakemake working directory with `--singularity-args "--bind <SNAKEMAKE_WORKING_DIRECTORY>"`

```sh
snakemake --sdm conda apptainer --singularity-args "--bind <SNAKEMAKE_WORKING_DIRECTORY>" --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going
```

### Run simulations in batches with singularity on slurm cluster

```sh
batches=50
for num in {1..$batches}
do
  snakemake --sdm conda apptainer --singularity-args "--bind <SNAKEMAKE_WORKING_DIRECTORY>" --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going --batch all_sims=$num/$batches all_sims
done
```

### Run workflow on local machine

The default snakemake profile is to run on a slurm cluster, but you can take any of the above commands and run snakemake on your local machine by adding `--profile profiles/local` to your snakemake command. Make sure to edit `workflow/profiles/local/config.yaml` to reflect the hardware limits of your local machine.

