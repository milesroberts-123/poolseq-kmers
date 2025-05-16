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
results = data.frame(
  jaccard = jaccard(kmerCounts$c.x, kmerCounts$c.y),
  bray_curtis = bray_curtis(kmerCounts$c.x, kmerCounts$c.y),
  cosine = cosine(kmerCounts$c.x, kmerCounts$c.y)
)

write.table(results, dissimOutput, row.names = F, quote = F, sep = "\t")
