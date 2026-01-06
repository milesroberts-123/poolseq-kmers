library("dplyr")

# number of replicates per simulation
replicates <- 1:3

# number of individuals to randomly sample at end of SLiM
# number of genomes in pool will be twice this number
sample_sizes <- c(25, 50, 75, 100, 125)

# sequencing coverage of entire pool = avg bp per read * num of reads / genome size
coverages <- c(10, 25, 50, 100, 150)

# insilico-seq error models
sequencers <- c("hiseq", "novaseq", "miseq", "nextseq")

# genome sizes
L <- c(1e6, 2e6, 4e6)

# k-mer sizes
k <- c(27, 31, 41)

# multinomial probabilities for generating k-mers
# must add to 1
pA <- 0.25
pT <- 0.25
pC <- 0.25
pG <- 0.25

# population parameters
population_sizes <- c(1000, 2000)

mutation_rates <- c(1e-8, 2e-8)

recombination_rates <- c(1e-8, 2e-8)

# shape parameter for zeta distribution
# large number = effectively no repeats
shapes <- c(2.854, 3.591, 4.772, 1e6)

# whether identical k-mer copies should be shuffled (TRUE) or concatenated (FALSE)
shuffles <- c(TRUE, FALSE)

# create data frame of workflow parameters
one_pop_params <- expand.grid(
  rep = replicates,
  N = population_sizes,
  n = sample_sizes,
  L = L,
  pA = pA,
  pC = pC,
  pG = pG,
  pT = pT,
  k = k,
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  sequencer = sequencers,
  simtype = "onepop",
  shape = shapes,
  shuffle = shuffles
)

# parameters for two population model
two_pop_params <- expand.grid(
  rep = replicates,
  N1 = population_sizes,
  N2 = population_sizes,
  mg1 = c(0),
  mg2 = c(0),
  tau = c(0.5, 1, 2),
  n = sample_sizes,
  L = L,
  pA = pA,
  pC = pC,
  pG = pG,
  pT = pT,
  k = k,
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  simtype = "twopop"
)

# convert tau in units of N generations to generations
two_pop_params$tau <- two_pop_params$tau * (two_pop_params$N1 + two_pop_params$N2)

sweep_params <- expand.grid(
  rep = replicates,
  N = population_sizes,
  n = sample_sizes,
  L = L,
  pA = pA,
  pC = pC,
  pG = pG,
  pT = pT,
  k = k,
  h = c(0, 0.5, 1),
  Nes = c(5, 10, 25, 50, 100),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  simtype = "sweep"
)

# convert Nes to selection coefficient
sweep_params$s <- sweep_params$Nes / sweep_params$N

params <- bind_rows(one_pop_params, sweep_params, two_pop_params)

# if needed, subset to one simulation type
params <- params[(params$simtype == "onepop"), ]

# replace all NA with 0, so that snakemake stays happy
params[is.na(params)] <- 0

# add simulation id
params$ID <- 1:nrow(params)

# add initial random seeds so that each simulation is completely reproducible
params$slimseed <- sample(1:((2^31) - 1), size = nrow(params), replace = F)
params$issseed <- sample(1:((2^31) - 1), size = nrow(params), replace = F)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
