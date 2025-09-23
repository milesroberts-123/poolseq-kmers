print("Loading libraries...")
library(VGAM)
library(ggplot2)

print("Parsing arguments...")
args = commandArgs(trailingOnly=TRUE)
pA = as.numeric(args[1])
pC = as.numeric(args[2])
pG = as.numeric(args[3])
pT = as.numeric(args[4])
shape = as.numeric(args[5])
k = as.numeric(args[6])
L = as.numeric(args[7])
ID = as.numeric(args[8])
#countKmers = as.numeric(args[8])
shuffleKmers = as.logical(args[9])

alphabet = c("A", "C", "G", "T")

print("Generating genome...")
custom_probs = c(pA, pC, pG, pT)

total_genome = c()

while(k*length(total_genome) < L){

  # generate random k-mer
  kmer = paste(sample(alphabet, size = k, replace = T, prob = custom_probs), collapse = "")

  # detmerine frequency of k-mer
  freq = rzeta(1, shape)

  # copy k-mer 
  kmers = rep(kmer, freq)

  # add k-mers to total_genome
  total_genome = c(total_genome, kmers)
}

# randomly order vector elements
if(shuffleKmers){
	print("Randomly ordering k-mers...")
	total_genome = sample(total_genome, replace = F)
}

# concatenate k-mers together
print("Concatenating shuffled k-mers...")
total_genome = paste(total_genome, collapse = "")

# trim genome to target length
print("Trimming genome to target length...")
total_genome = substring(total_genome, 1, L)

# write genome output
#print(total_genome)
print("Writing genome to fasta file...")
writeLines(c("> 1", total_genome), paste("ancestral_genome_results/", ID, ".fasta", sep = ""))

# count k-mers in synthetic genome to confirm zipf distribution
#print("Counting k-mers genome...")
#start = 1
#stop = k
#
#kmer_list = list()
#
#while(stop <= L){
#  
#  kmer = substring(total_genome, start, stop)
#  
#  if(is.null(kmer_list[[kmer]])){
#      kmer_list[[kmer]] = 1
#  } else {
#    kmer_list[[kmer]] = kmer_list[[kmer]] + 1
#  }
#
#  start = start + 1
#  stop = stop + 1
#}
#
## plot zipf distribution
#print("Plotting power-law distribution...")
#kmer_freqs = table(unlist(kmer_list))
#
#plotdata = data.frame(
#  freqs = as.numeric(names(kmer_freqs)),
#  ranks = kmer_freqs
#)
#
#ggplot(plotdata, aes(x = log10(freqs), y = log10(ranks.Freq))) +
#  geom_point() +
#  theme_classic() +
#  labs(x = "k-mer count", y = "Number of unique k-mers with count")
#  
#ggsave(paste("power_law_", ID, ".jpg", sep = ""), dpi = 350)
