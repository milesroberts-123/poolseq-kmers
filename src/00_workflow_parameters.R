# create data frame of workflow parameters
params = data.frame(
  N = 1000,
  sigma = 0,
  n = 25,
  mu = 1e-8,
  R = 1e-8,
  kappa = 10000,
  cov = 100,
  L = 1e6
)

# save
write.table(params, "../config/paramters.tsv", sep = "\t", quotes = F, row.names = F)