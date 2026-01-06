# Resources

This directory contains two resources:

1. An example R script for generating config/parameters.tsv

2. An example slurm script to launch the snakemake workflow on a slurm cluster. You submit this script with `sbatch` and the script will run an instance of snakemake, which will itself go on to submit and manage additional slurm jobs for each step in the workflow. You will want to change the SBATCH header lines depending on the configuration of your specific HPCC system. The script also includes examples of (a) running the whole snakemake workflow at once with conda, (b) running the snakemake workflow in batches with conda, (c) running the whole snakemake workflow with singularity, (d) running the snakemake workflow in batches with singularity
