# create data frame of workflow parameters
params = expand.grid(
  rep = c(1:10),
  N = c(1000),
  n = c(5, 10, 25, 50, 100),
  sigma = c(0),
  mu = c(1e-8),
  R = c(1e-8),
  cov = c(20, 50, 80, 100, 200),
  L = c(1e6)
)

params$ID = 1:nrow(params)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
