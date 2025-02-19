library("dplyr") 

replicates = 1:5
sample_sizes = c(40, 50, 75, 100, 125)
coverages = c(50, 100, 200, 250, 300)
sequencers = c("miseq", "hiseq", "nextseq", "novaseq")

#sample_sizes = c(5, 10)
#coverages = c(20, 50)
#sequencers = c("miseq", "hiseq")

# create data frame of workflow parameters
params = expand.grid(
  rep = replicates,
  N = c(1000),
  n = sample_sizes,
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = coverages,
  L = c(1e6),
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
#two_pop_params = expand.grid(
#  rep = replicates,
#  N1 = c(1000),
#  N2 = c(1000),
#  mg1 = c(0, 0.1, 0.2),
#  mg2 = c(0, 0.1, 0.2),
#  n = sample_sizes,
#  sigma = c(0),
#  mu = c(1e-8),
#  R = c(1e-8),
#  cov = coverages,
#  L = c(1e6),
#  sequencer = sequencers,
#  simtype = "twopop"
#)

# selective sweep parameters
#sweep_params = expand.grid(
#  rep = replicates,
#  N = c(1000),
#  n = sample_sizes,
#  h = c(0, 0.5, 1),
#  s = c(0.1, 0.25, 0.5),
#  sigma = c(0),
#  mu = c(1e-8),
#  R = c(1e-8),
#  cov = coverages,
#  L = c(1e6),
#  sequencer = sequencers,
#  simtype = "sweep"
#)

# combine all parameters into one table
#params = bind_rows(one_pop_params, two_pop_params, sweep_params)

# add simulation id
params$ID = 1:nrow(params)

# replace all NA with 0, so that snakemake stays happy
params[is.na(params)] = 0

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
