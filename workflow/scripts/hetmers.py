# I am using the algorithm published in smudgeplot
# https://github.com/KamilSJaron/smudgeplot/blob/master/playground/more_away_pairs.py
#from Bio.Seq import Seq
import pandas as pd
import numpy as np
import hashlib
import click
import datetime

@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", default=None, help="Path to input k-mer count table (KMC or Jellyfish)", multiple=False)
@click.option("-a", "--alleles", default=[2], help="Number of alleles", multiple=True)
@click.option("-m", "--minimum", default=5, help="Minimum required k-mer count", multiple=False)
@click.option("-o", "--output-prefix", default="hetmers", help="Prefix for output files")
def main(input, alleles, minimum, output_prefix):
  """
  Script similar to smudgeplot hetkmers, except we keep the hetmer sequences
  """

  print(datetime.datetime.now(), ": Loading k-count file...")
  counts = pd.read_csv(input, sep="\t", names = ["seq", "count"])

  print(datetime.datetime.now(), ": Filtering low-frequency k-mers...")
  filtered_counts = counts[counts['count'] >= minimum]
  del counts

  # separate seqs and counts
  seqs = filtered_counts['seq'].values
  counts = filtered_counts['count'].values
  del filtered_counts

  # determine k based on length of first string
  k = len(seqs[1])
  print(datetime.datetime.now(), ": K is " + str(k))

  print(datetime.datetime.now(), ": Extracting all but middle base...")
  k_half = k//2
  no_center_ks = [s[:k_half]+s[-k_half:] for s in seqs]

  # ger reverse complement
  print(datetime.datetime.now(), ": Reverse complementing...")
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  no_center_ks_rev = ["".join(complement.get(base, base) for base in reversed(x)) for x in no_center_ks]

  # has both forward and reverse complement
  print(datetime.datetime.now(), ": Hashing...")
  no_center_ks_hash = [int(hashlib.sha256(x.encode('utf-8')).hexdigest(), 16) for x in no_center_ks]
  no_center_ks_rev_hash = [int(hashlib.sha256(x.encode('utf-8')).hexdigest(), 16) for x in no_center_ks_rev]
  del no_center_ks
  del no_center_ks_rev

  # use the smaller of the two hashes as the sequence id
  print(datetime.datetime.now(), ": Getting the minimum hash...")
  min_hashs = [min(x,y) for x,y in zip(no_center_ks_hash,no_center_ks_rev_hash)]
  del no_center_ks_hash
  del no_center_ks_rev_hash

  # group kmers based on hash
  print(datetime.datetime.now(), ": Grouping unique hashes into a dictionary...")
  d = {}
  for i, num in enumerate(min_hashs):
      if num in d:
          d[num].append(i)
      else:
          d[num] = [i]

  # Filter to keep only the numbers with more than one occurrence
  print(datetime.datetime.now(), ": Filter hashes...")
  ans = {key: value for key, value in d.items() if len(value) in alleles}
  del d

  # extract het mer sequences and counts
  print(datetime.datetime.now(), ": Extracting counts and sequences...")
  hetmers = [seqs[value] for key, value in ans.items()]
  hetmer_counts = [counts[value] for key, value in ans.items()]

  # save results
  print(datetime.datetime.now(), ": Saving results...")
  np.savetxt(output_prefix + "_seqs.csv", hetmers, delimiter = ",", fmt='%s')
  np.savetxt(output_prefix + "_counts.csv", hetmer_counts, delimiter = ",")

  print(datetime.datetime.now(), ": Done! :D")

if __name__ == '__main__':
    main()


