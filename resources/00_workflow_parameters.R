library("dplyr") 

replicates = 1:3
sample_sizes = c(25, 50, 75, 100)
coverages = c(50, 100, 150, 200)
sequencers = c("miseq")
#sequencers = c("miseq", "hiseq", "nextseq", "novaseq")

population_sizes = c(1000)
mutation_rates = c(1e-8)
recombination_rates = c(1e-8)

shapes = c(4.22, 1e6)
shuffles = c(TRUE, FALSE)
chrom_length = 2e6

# create data frame of workflow parameters
one_pop_params = expand.grid(
  rep = replicates,
  #rep = replicates,
  N = population_sizes,
  n = sample_sizes,
  sigma = c(0),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  L = chrom_length,
  sequencer = sequencers,
  simtype = "onepop",
  shape = shapes,
  shuffle = shuffles,
  pA = 0.25,
  pC = 0.25,
  pG = 0.25,
  pT = 0.25
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
  L = chrom_length,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  pA = 0.25,
  pC = 0.25,
  pG = 0.25,
  pT = 0.25,
  simtype = "twopop"
)

# convert tau in units of N generations to generations
two_pop_params$tau = two_pop_params$tau*(two_pop_params$N1 + two_pop_params$N2)

# selective sweep parameters
sweep_params = expand.grid(
  rep = replicates,
  #rep = replicates,
  N = population_sizes,
  n = sample_sizes,
  h = c(0, 0.5, 1)
  #h = c(0.5),
  Nes = c(10, 25, 50, 100),
  #s = 0.1,
  sigma = c(0),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  L = chrom_length,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  pA = 0.25,
  pC = 0.25,
  pG = 0.25,
  pT = 0.25,
  simtype = "sweep"
)

# convert Nes to selection coefficient
sweep_params$s = sweep_params$Nes/sweep_params$N

# bulk-segregant analysis parameters
bsa_params = expand.grid(
  rep = replicates,
  N = 2000,
  n = sample_sizes,
  sigma = c(0),
  mu = mutation_rates,
  R = recombination_rates,
  cov = coverages,
  L = chrom_length,
  sequencer = sequencers,
  shape = shapes,
  shuffle = shuffles,
  pA = 0.25,
  pC = 0.25,
  pG = 0.25,
  pT = 0.25,
  qtl_mean = 0,
  qtl_sigma = 1,
  qtl_prop = 0.1,
  optimum_mean = 0,
  optimum_sigma = 10,
  phenotype_cutoff = 0.05,
  simtype = "bsa"
)

# combine all parameters into one table
params = bind_rows(one_pop_params, two_pop_params, sweep_params, bsa_params)
#params = bind_rows(one_pop_params, two_pop_params)
#params = bind_rows(one_pop_params, sweep_params)
#params = bind_rows(sweep_params, bsa_params)
#params = sweep_params

# replace all NA with 0, so that snakemake stays happy
params[is.na(params)] = 0

# subset if needed
params = params[(params$simtype %in% c("sweep", "bsa")),]

# add simulation id
params$ID = 1:nrow(params)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)

# generate deletions in reference genome
#deletions = rgeom(120000, 0.2)
#
#deletions = deletions[(deletions > 0)]
#
#ends = cumsum(deletions)
#starts = c(0, head(ends,-1))
#
#mybed = data.frame(
# chrom = 1,
# starts = starts,
# ends = ends
#)
#
#write.table(mybed, "../config/mask.bed", sep = "\t", quote = F, row.names = F, col.names = F)
