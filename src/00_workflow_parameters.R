library("dplyr") 

replicates = 1
#sample_sizes = c(25, 50, 75, 100, 125, 150)
#coverages = c(50, 100, 150, 200, 250, 300)
#sequencers = c("miseq", "hiseq", "nextseq", "novaseq")

sample_sizes = c(50, 75, 100, 125)
coverages = c(50, 100, 150, 200)
sequencers = c("miseq", "hiseq", "nextseq", "novaseq")

chrom_length = 2e6

# create data frame of workflow parameters
one_pop_params = expand.grid(
  rep = replicates,
  N = c(1000),
  n = sample_sizes,
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = chrom_length,
  sequencer = sequencers,
  simtype = "onepop",
  N1 = 0,
  N2 = 0,
  mg1 = 0,
  mg2 = 0,
  h = 0,
  s = 0
)

# parameters for two population model
two_pop_params = expand.grid(
  rep = replicates,
  N1 = c(1000),
  N2 = c(1000),
  mg1 = c(0),
  mg2 = c(0),
  tau = c(2000),
  n = sample_sizes,
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = chrom_length,
  sequencer = sequencers,
  simtype = "twopop"
)

# selective sweep parameters
sweep_params = expand.grid(
  rep = replicates,
  N = c(1000),
  n = sample_sizes,
  h = 0.5,
  s = 0.5,
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = chrom_length,
  sequencer = sequencers,
  simtype = "sweep"
)

# combine all parameters into one table
params = bind_rows(one_pop_params, two_pop_params, sweep_params)
#params = bind_rows(one_pop_params, two_pop_params)

# add simulation id
params$ID = 1:nrow(params)

# replace all NA with 0, so that snakemake stays happy
params[is.na(params)] = 0

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)

# generate deletions in reference genome
deletions = rgeom(50000, 0.1)

deletions = deletions[(deletions > 0)]

ends = cumsum(deletions)
starts = c(0, head(ends,-1))

mybed = data.frame(
 chrom = 1,
 starts = starts,
 ends = ends
)

write.table(mybed, "../config/mask.bed", sep = "\t", quote = F, row.names = F, col.names = F)
