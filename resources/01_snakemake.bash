#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --partition=savio2_htc
#SBATCH --account=co_moilab

# output information about how this job is running using bash commands
echo "This job is running on $HOSTNAME on `date`"

# Load conda module, helps nodes find my conda path for some reason
module purge
module load Conda/3

# load snakemake
echo Loading snakemake...
conda activate snakemake-NEW

# go to workflow directory with Snakefile
echo Changing directory...
cd ../workflow

# unlock snakemake if previous instance of snakemake failed
echo Unlocking snakemake...
snakemake --unlock --cores 1 --batch all=1/100

# submit snakemake to HPCC
# subtract one job and one core from max to account for this submission command
# rerun-incomplete in case previous snakemake instances failed and left incomplete files
# Max cpu count for my SLURM account is 1040, subtract 1 to account for scheduler
# Max job submit count is 1000, subtract 1 to account for scheduler

## RUN WHOLE WORKFLOW ON SLURM CLUSTER WITH CONDA ##

#snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going

## RUN WHOLE WORKFLOW ON SLURM CLUSTER WITH SINGULARITY + CONDA ##

#snakemake --sdm conda apptainer --singularity-args "--bind /mnt/scratch/robe1195/Josephs_Lab_Projects/poolseq-kmers/workflow" --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going

## RUN WORKFLOW IN BATCHES ON SLURM CLUSTER WITH CONDA ##

#for num in {1..50}
#do
#  snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going --batch all=$num/50
#done

## RUN WORKFLOW IN BATCHES ON SLURM CLUSTER WITH SINGULARITY + CONDA ##

batch=10
for i in $( eval echo {1..$batch} )
do
  snakemake --sdm conda apptainer --singularity-args "--bind /global/scratch/users/milesroberts/poolseq-kmers/workflow/" --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going --batch all=$i/$batch
done
