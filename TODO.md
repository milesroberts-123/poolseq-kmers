# To do

## high priority

- [ ] run final batch of simulations

- [ ] add bcftools call as another option after vg

- [x] run histogram branch for final time

- [x] add varscan to empirical data workflow

- [x] polish freqk, skip over variants close to chromosome ends, non-isolated variants

- [x] add workflow for processing empirical datasets

- [x] polish simulation repo: remove left overs

- [x] polish freqk: separate code into modules, skip over problematic variants (close to chromosome ends)

- [x] have vg in simulations also map unpaired reads

- [x] add workflow to compare freqk and vg on empirical data from plantpan

- [x] analyze zeta distribution across ncbi reference genomes

- [x] use zeta distribution to inform parameter choices for ancestral genome generation

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

- [x] change wildcards to something like: {simulation id}_{population id} so that I don't need separate rules for 1 population vs 2 population workflows?

- [x] [add minimum snakemake version](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#depend-on-a-minimum-snakemake-version)

- [x] freebayes

- [x] add seeds to iss and slim so that unit tests will always give same answer

- [x] add dissimilarity script back in. I couldn't figure this out - even when using the branch function.

- [x] [add snakefmt via github actions](https://github.com/snakemake/snakefmt?tab=readme-ov-file#github-actions)

- [x] add freqk

## lower priority

- [ ] add structural variants randomly to samples VCF file output from slim

- [ ] add angsd?

- [ ] add time series slim simulation

- [ ] add `bcftools call`?

- [ ] calculate dxy from slim outputs

- [ ] add job groups?

- [ ] add more sequencing simulators: dwgsim, mason, or add another sequencer error profile

- [ ] [Try rewriting discosnp as a shadow rule](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#shadow-rules)

- [ ] add R notebook to snakemake

- [ ] write hetmers binary to calculate fst

- [ ] unit tests

- [ ] integration tests

- [ ] github actions

- [ ] generalize bayes theorem to negative binomial distribution

- [ ] snakemake reports

- [ ] add in ploidyfrost

- [ ] try adding kmer2snp

- [ ] add purifying selection simulation - does this also create unitigs?

- [ ] update to smudgeplot >0.3.0, once we're able to get k-mer sequences again

- [ ] figure out how to use unpaired reads in varscan

