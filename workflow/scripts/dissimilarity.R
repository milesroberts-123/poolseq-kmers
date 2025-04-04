rm(list = ls())

# get file names passed to script as arguments
print("Parsing arguments...")
args = commandArgs(trailingOnly=TRUE)
kmerFile1 = args[1]
kmerFile2 = args[2]
dissimOutput = args[3]

print("Input files:")
print(kmerFile1)
print(kmerFile2)

print("Output files:")
print(dissimOutput)

# dissimilarity functions
jaccard = function(x,y){
  union = length(x)
  
  inter = sum( (x > 0) & (y > 0) )
  
  1 - inter/union
}

bray_curtis = function(x,y){
  x = x/sum(x)
  y = y/sum(y)
  1 - 2*sum(pmin(x,y))/sum(x + y)
}

cosine = function(x,y){
  num = sum(x*y)
 
  a = sqrt(sum(x^2))
  b = sqrt(sum(y^2))

  1 - num/(a*b)
}

kmer_fst = function(x,y){
  x = x/sum(x)
  y = y/sum(y)
  
  xsq = x^2
  ysq = y^2
  
  pi_wn = 0.5 * ( (1-sum(xsq)) + (1-sum(ysq)) )
  
  pi_bw = 1-sum(x*y)
  
  pi_to = 0.5*pi_wn + 0.5*pi_bw
  
  fst_nei = 1 - pi_wn/pi_to
  
  fst_hud = 1- pi_wn/pi_bw
  
  return(c(fst_nei, fst_hud))
}

print("Reading files into memory...")
kmerCounts1 = read.table(kmerFile1, header = T, col.names = c("k", "c"), sep = "\t")
kmerCounts2 = read.table(kmerFile2, header = T, col.names = c("k", "c"), sep = "\t")

head(kmerCounts1)
head(kmerCounts2)

print("Merging files by k-mer...")
kmerCounts = merge(kmerCounts1, kmerCounts2, by = "k", all.x = T, all.y = T)

head(kmerCounts)

# change NAs to 0
print("Changing NAs to 0s...")
kmerCounts[is.na(kmerCounts)] = 0

# calculate k-mer dissimilarity metrics
print("Calculating dissimilarity...")

fst_values = kmer_fst(kmerCounts$c.x, kmerCounts$c.y)

results = data.frame(
  jaccard = jaccard(kmerCounts$c.x, kmerCounts$c.y),
  bray_curtis = bray_curtis(kmerCounts$c.x, kmerCounts$c.y),
  cosine = cosine(kmerCounts$c.x, kmerCounts$c.y),
  nei_fst = fst_values[1],
  hudson_fst = fst_values[2]
)

write.table(results, dissimOutput, row.names = F, quote = F, sep = "\t")
