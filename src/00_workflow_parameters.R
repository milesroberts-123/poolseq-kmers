# create data frame of workflow parameters
params = data.frame(
  ID = 1:5,
  N = 1000,
  sigma = 0,
  n = c(5, 10, 25, 50, 100),
  mu = 1e-8,
  R = 1e-8,
  kappa = 10000,
  cov = 100,
  L = 1e6
)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
