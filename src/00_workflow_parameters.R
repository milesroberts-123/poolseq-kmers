# create data frame of workflow parameters
params = expand.grid(
  rep = c(1:5),
  N = c(500, 1000, 5000, 10000),
  n = c(5, 10, 25, 50, 100),
  sigma = c(0, 0.5, 1),
  mu = c(1e-9, 5e-9, 1e-8),
  R = c(1e-9, 1e-8, 1e-7),
  cov = c(20, 50, 80, 100, 200),
  L = c(5e5, 1e6, 5e6)
)

params$ID = 1:nrow(params)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
