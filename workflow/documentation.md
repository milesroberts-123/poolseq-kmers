
# Overview

In brief, for each simulation this workflow does the following tasks:

* Create a ancestral genome 

* Run a slim simulation on the ancestral genome

* Simulate whole genome sequencing of a pool of genomes 

* Perform basic read quality control on the reads

* count k-mers and call SNPs in the reads

* Output read count and het-mer count statistics
 
See `images/rulegraph.svg` for a full example of how rules connect together for each simulation.

# Rules

## ancestral_genome.smk

## bcftools_discosnp.smk

## bcftools_fill_tags.smk

## bcftools_freebayes.smk

## bcftools_get_samples.smk

## bcftools_poolsnp.smk

## bcftools_remove_ref.smk

## blast.smk

## bwa_index.smk

## bwa_mem.smk

## common.smk

## compress.smk

## discosnp.smk

## dissimilarity.smk

## fastp.smk

## freebayes.smk

## hetmers.smk

## iss.smk

## kmc.smk

## makeblastdb.smk

## multiqc.smk

## pankmer.smk

## poolsnp.smk

## samtools_faidx.smk

## seqkit_get_ref.smk

## seqkit_get_samples.smk

## seqkit_hetmer_key.smk

## seqkit_rename.smk

## slim.smk
## smudgeplot.smk
## unitig_caller.smk
## varscan.smk
