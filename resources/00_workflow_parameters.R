library("dplyr") 

replicates = 1:5
#replicates = 1
#sample_sizes = c(25, 50, 75, 100, 125, 150)
sample_sizes = c(25, 50, 75)
#coverages = c(30, 50, 100, 150, 200, 250)
coverages = c(50, 100, 150)
#sequencers = c("hiseq", "novaseq", "miseq", "nextseq")
sequencers = "hiseq"

population_sizes = c(1000)
mutation_rates = c(1e-8)
recombination_rates = c(1e-8)

shapes = c(2.5, 5, 10, 15, 1e6)
shuffles = c(TRUE, FALSE)

# create data frame of workflow parameters
one_pop_params = expand.grid(
  rep = replicates,
  N = population_sizes,
  n = sample_sizes,
  sigma = c(0),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  sequencer = sequencers,
  simtype = "onepop",
  shape = shapes,
  shuffle = shuffles
)

# parameters for two population model
two_pop_params = expand.grid(
  rep = replicates,
  N1 = population_sizes,
  N2 = population_sizes,
  mg1 = c(0),
  mg2 = c(0),
  tau = c(0.5, 1, 2),
  n = sample_sizes,
  sigma = c(0),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  simtype = "twopop"
)

# convert tau in units of N generations to generations
two_pop_params$tau = two_pop_params$tau*(two_pop_params$N1 + two_pop_params$N2)

sweep_params = expand.grid(
  rep = replicates,
  N = population_sizes,
  n = sample_sizes,
  #h = c(0, 0.5, 1),
  h = 0.5,
  Nes = c(5, 10, 25, 50, 100),
  sigma = c(0),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  simtype = "sweep"
)

# convert Nes to selection coefficient
sweep_params$s = sweep_params$Nes/sweep_params$N

params = bind_rows(one_pop_params, sweep_params, two_pop_params)

# if needed, subset to one simulation type
params = params[(params$simtype == "sweep"),]

# replace all NA with 0, so that snakemake stays happy
params[is.na(params)] = 0

# add simulation id
params$ID = 1:nrow(params)

# add initial random seeds
params$slimseed = sample(1:((2^31)-1), size = nrow(params), replace = F)
params$issseed = sample(1:((2^31)-1), size = nrow(params), replace = F)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
