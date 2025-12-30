# poolseq-kmers

[![Super-Linter](https://github.com/milesroberts-123/poolseq-kmers/actions/workflows/linter.yml/badge.svg)](https://github.com/marketplace/actions/super-linter)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>)

Author: Miles Roberts

## Table of Contents

[Overview](#overview)

[Setup](#setup)

[Inputs](#inputs)

[Outputs](#outputs)

[Running the workflow](#running-the-workflow)

[Statistical analysis and figure creation](#statistical-analysis-and-figure-creation)

[Citations](#citations)
 
## Overview

Simulation workflow to investigate the utility of k-mers, het-mers, and k-unitigs for pool-seq data analysis built with snakemake (v 9.3.3) and run with the snakemake slurm plugin (v 1.3.6) and snakedeploy (v 0.11.0)

## Setup

1. Install conda/mamba

2. Download the workflow from github

3. Create a mamba environment with snakemake and any plugins you need

This is how to make new mamba environment named snakemake with snakemake and the slurm plugin installed. If you are not running the workflow on a SLURM cluster, you can install a different pluggin

```
conda create -y -f snakemake-mamba-env.yaml

conda activate snakemake
```

4. Check the snakemake profile for the proper executer. The default profile runs snakemake on a slurm cluster (`workflow/profiles/default/config.yaml`), but you should still change the slurm account, slurm partition, and default resources to match your system.

## Inputs

See config/README.md for a complete description of workflow inputs. In short, you need two files:

* config/config.yaml: describes parameters that are held constant for every job in the workflow

* config/parameters.tsv: is a table of parameters that can vary between simulations. Each simulation corresponds to a different row, and each parameter is a column. Each simulation should have a column `ID` that is an integer used as a unique identifier.

## Outputs

For each simulation, this workflow outputs:

* SNP calls from Varscan

* SNP calls from PoolSNP

* SNP calls from discosnp

* Het-mers from smudgeplot

* Het-mers from hetmers

## Usage

Examples commands are in `resources/01_snakemake.bash`

### Run whole workflow with conda envs on slurm cluster

`snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going`

### Run workflow in batches with conda envs on slurm cluster

```
for num in {1..50}
do
  snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going --batch all=$num/50
done
```

*Tip:* Collapse logs folder into one archive to minimize number of files on your system

```
# create initial archive
tar -cvf logs.tar logs/

# add more files
tar -uvf logs.tar logs/

# compress at the very end
gzip logs.tar
```

### Run workflow whole workflow at once with singularity on slurm cluster

Instead of downloading and building all of the conda environments, you can just download a container with all of the conda environments pre-installed.

Need to pass `--use-singularity` to snakemake and also your snakemake working directory with `--singularity-args "--bind <SNAKEMAKE_WORKING_DIRECTORY>"`

```
snakemake --sdm conda apptainer --singularity-args "--bind ~/Josephs_Lab_Projects/poolseq-kmers/workflow" --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going
```

### Run workflow in batches with singularity on slurm cluster

```
for num in {1..50}
do
  snakemake --sdm conda apptainer --singularity-args "--bind /mnt/scratch/robe1195/Josephs_Lab_Projects/poolseq-kmers/workflow" --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going --batch all=$num/50
done
```

### Run workflow on local machine

The default snakemake profile is to run on a slurm cluster, but you can take any of the above commands and run snakemake on your local machine by adding `--profile profiles/local` to your snakemake command. Make sure to edit `workflow/profiles/local/config.yaml` to reflect the hardware limits of your local machine.

## Citations

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

### Adding rule-specific resources to profile

https://github.com/snakemake/snakemake-executor-plugin-slurm/blob/main/docs/further.md

### Visualize DAG

`snakemake --dag | dot -Tpdf > dag.pdf`

### Snakemake reports

Could be a good idea to move my notebook into scripts, which will generate figures that get compiled into a snakemake report. In my mind, this is more easily reproducibile than having someone configure a notebook. However, you can't really explore the data this way beyond whatever figures you predetermine.

Another option is to just add my notebook as a snakemake rule:

https://nbis-reproducible-research.readthedocs.io/en/course_2104/rmarkdown/#r-markdown-and-snakemake

This could be nicer because R markdown will give me more control over what the report.html looks like.
