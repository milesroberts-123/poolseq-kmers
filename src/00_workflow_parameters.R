library("dplyr") 

replicates = 1
sample_sizes = c(5, 10, 25, 50, 100)
coverages = c(20, 50, 80, 100, 200)
sequencers = c("miseq", "hiseq", "nextseq", "novaseq")

# create data frame of workflow parameters
one_pop_params = expand.grid(
  rep = replicates,
  N = c(1000),
  n = sample_sizes,
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = c(1e6),
  sequencer = sequencers,
  simtype = "onepop"
)

# parameters for two population model
two_pop_params = expand.grid(
  rep = replicates,
  N1 = c(1000),
  N2 = c(1000),
  mg1 = c(0, 0.1, 0.2, 0.3),
  mg2 = c(0, 0.1, 0.2, 0.3),
  n = sample_sizes,
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = c(1e6),
  sequencer = sequencers,
  simtype = "twopop"
)

# selective sweep parameters
sweep_params = expand.grid(
  rep = replicates,
  N = c(1000),
  n = sample_sizes,
  h = c(0.1, 0.5, 0.9),
  s = c(0.1, 0.2, 0.3),
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = c(1e6),
  sequencer = sequencers,
  simtype = "sweep"
)

# combine all parameters into one table
params = bind_rows(one_pop_params, two_pop_params, sweep_params)

# add simulation id
params$ID = 1:nrow(params)

# replace all NA with 0, so that snakemake stays happy
params[is.na(params)] = 0

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
