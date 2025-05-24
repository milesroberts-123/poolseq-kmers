# poolseq-kmers

Author: Miles Roberts

Simulation workflow to investigate the utility of k-mers, het-mers, and k-unitigs for pool-seq data analysis built with snakemake (v 9.3.3) and run with the snakemake slurm plugin (v 1.3.6) and snakedeploy (v 0.11.0)

## Table of Contents

* Overview

* Setup

* Inputs

* Outputs

* Running the workflow

* Statistical analysis and figure creation

## Overview

## Setup

1. Install mamba

2. Create a mamba environment with snakemake and any plugins you need

This is how to make new mamba environment named snakemake with snakemake and the slurm plugin installed. If you are not running the workflow on a SLURM cluster, you can install a different pluggin

```
mamba create -y -n snakemake snakemake snakemake-executor-plugin-slurm snakedeploy

mamba activate snakemake
```

3. Download the workflow from github

4. Check the snakemake profile for the proper executer. The default profile runs snakemake on a slurm cluster (`workflow/profiles/default/config.yaml`), but you should still change the slurm account, slurm partition, and default resources to match your system.

## Inputs

See config/README.md for a complete description of workflow inputs. In short, you need two files:

* config/config.yaml: describes parameters that are held constant for every simulation in the workflow

* config/parameters.tsv: is a table of parameters that vary between simulations. Each simulation corresponds to a different row, and each parameter is a column. Each simulation should have a column `ID` that is an integer used as a unique identifier.

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

`snakemake --cluster "sbatch --time={resources.time} --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb_per_cpu} --partition=josephsnodes --account=josephsnodes --output=logs/slurm/%j.out --error=logs/slurm/%j.out" --jobs 950 --cores 950 --use-conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going`

### Run workflow in batches with conda envs on slurm cluster

```
for num in {1..50}
do
  snakemake --cluster "sbatch --time={resources.time} --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb_per_cpu} --partition=josephsnodes --account=josephsnodes --output=logs/slurm/%j.out --error=logs/slurm/%j.out" --jobs 975 --cores 975 --use-conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going --batch all=$num/50
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
snakemake --sdm conda apptainer --singularity-args "--bind ~/Josephs_Lab_Projects/poolseq-kmers/workflow" --cores 1
```

### Run workflow in batches with singularity on slurm cluster

### Run workflow on local machine

`snakemake --profile profiles/local`

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

- [x] upgrade to latest snakemake version

- [x] add slurm profile

- [x] add rule-specific resources to profile

- [x] [Add workflow hub requirements](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#uploading-workflows-to-workflowhub)

- [x] add local profile

- [x] try calculating fst in R with hetmers from individual pools and combined pools -> this will pilot my idea before I try coding it into rust

- [x] use lookup functions

- [x] split seqkit rules into two rules

- [x] snakefmt

- [x] add local rules

- [x] parameter space coding

- [ ] [Try rewriting discosnp as a shadow rule](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#shadow-rules)

- [ ] add R notebook to snakemake

- [ ] unit tests

- [ ] integration tests

- [ ] github actions

- [ ] write hetmers binary to calculate fst

### lower priority

- [ ] generalize bayes theorem to negative binomial distribution

- [ ] snakemake reports

- [ ] add in ploidyfrost

- [ ] try adding kmer2snp

- [ ] add purifying selection simulation - does this also create unitigs?

- [ ] update to smudgeplot >0.3.0, once we're able to get k-mer sequences again

- [ ] figure out how to use unpaired reads in varscan 
