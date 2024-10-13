# I am using the algorithm published in smudgeplot
# https://github.com/KamilSJaron/smudgeplot/blob/master/playground/more_away_pairs.py

# define a hash function such that a k-mer and it's reverse complement return the same value
# https://bioinformatics.stackexchange.com/questions/3486/how-to-write-a-hash-function-for-canonical-kmers
def kmer_to_int:

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
