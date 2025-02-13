# I am using the algorithm published in smudgeplot
# https://github.com/KamilSJaron/smudgeplot/blob/master/playground/more_away_pairs.py
from Bio.Seq import Seq

# parameters
# i: input
# o: output
# c: k-mer count threshold
# m: 2,3,4 whether to consider bi-allelic, tri-allelic, or tetra-allelic sites

# define a hash function such that a k-mer and it's reverse complement return the same value
# https://bioinformatics.stackexchange.com/questions/3486/how-to-write-a-hash-function-for-canonical-kmers
def kmer_hash(x):
  # get reverse complement
  x_rev = x.reverse_complement()
  
  # convert to strings
  x_str = str(x)
  x_rev_str = str(x_rev)
  
  # hash strings
  for_hash = hash(x_str)
  rev_hash = hash(x_rev_str)
  
  # return minimum hash
  return(min(for_hash, rev_hash))
  
def main:
  # read k-mer count file

  # split k-mers into halves, except the center bp

  # hash k-mer halves

  # get list of all unique hash values

  # loop over each unique hash value, group k-mers according to hash value
  pairs = {}
  for hash_value in hash_values:
    kmerfile[kmerfile["hash"] == value]
    
  # calculate allele frequencies
  
  # for two populations, calculate Fst





# split k-mer in half
k_l = k//2
k_r = k - k_l

#initialize dictionaries in which the key is the hash of half of the kmer, and the value is a list of indices of the kmers with that same hash
kmer_L_hashes = defaultdict(list)
kmer_R_hashes = defaultdict(list)

#initialize pairs, which will be returned by get_1away_pairs
pairs = []

#initialize dictionaries containing the left halves and the right halves (since we will have to check cases where the left half differs by 1 and the right half differs by 1)
local_index_to_kmer_L = {}
local_index_to_kmer_R = {}

#for each kmer, calculate its left hash and right hash, then add its index to the corresponding entries of the dictionary
for i, kmer in local_index_to_kmer.items():
  kmer_L = kmer[:k_L]
  kmer_R = kmer[k_L:]
  local_index_to_kmer_L[i] = kmer_L
  local_index_to_kmer_R[i] = kmer_R
  kmer_L_hashes[kmer_to_int(kmer_L)] += [i]
  kmer_R_hashes[kmer_to_int(kmer_R)] += [i]

#for each left hash in which there are multiple kmers with that left hash, find the list of pairs in which the right half differs by 1. (aka, if left half matches, recurse on right half).
for kmer_L_hash_indices in kmer_L_hashes.values(): #same in first half
  if len(kmer_L_hash_indices) > 1:
    pairs += 

#for each right hash in which there are multiple kmers with that right hash, find the list of pairs in which the left half differs by 1. (aka, if right half matches, recurse on left half).
for kmer_R_hash_indices in kmer_R_hashes.values(): #same in second half
  if len(kmer_R_hash_indices) > 1:
    pairs += 
