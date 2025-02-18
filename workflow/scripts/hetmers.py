# I am using the algorithm published in smudgeplot
# https://github.com/KamilSJaron/smudgeplot/blob/master/playground/more_away_pairs.py
#from Bio.Seq import Seq
import pandas as pd
import numpy as np
import hashlib
import click
import datetime

# function to get hetmers
def get_hetmers(input, type, alleles, minimum, output_prefix):
  """
  Script similar to smudgeplot hetkmers, except we retain the hetmer sequences.
  """
  print(datetime.datetime.now(), ": Loading k-mer count file " + input + "...")
  counts = pd.read_csv(input, sep="\t", names = ["seq", "count"], header = None)

  print(datetime.datetime.now(), ": Filtering k-mers with count less than " + minimum + "...")
  filtered_counts = counts[counts['count'] >= minimum]
  del counts

  # separate seqs and counts
  seqs = filtered_counts['seq'].values
  counts = filtered_counts['count'].values
  del filtered_counts

  # determine k based on length of first string
  k = len(seqs[1])
  print(datetime.datetime.now(), ": k is " + str(k))

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
  del seqs
  del counts
  del ans

  # save results
  print(datetime.datetime.now(), ": Saving results...")
  np.savetxt(output_prefix + "_seqs.csv", hetmers, delimiter = ",", fmt='%s')
  np.savetxt(output_prefix + "_counts.csv", hetmer_counts, delimiter = ",")

  print(datetime.datetime.now(), ": Done! :D")

# function to work with one hetmer table
def get_allele_freqs(hetmer_table, output_prefix):
  print(datetime.datetime.now(), ": Loading hetmer counts from " + hetmer_table + "_counts.csv" +  "...")
  counts = pd.read_csv(hetmer_table + "_counts.csv", header=None)

  print(datetime.datetime.now(), ": Calculating frequencies...")
  a = counts.iloc[:,0].to_numpy()
  A = counts.iloc[:,1].to_numpy()
  freqs = a/(a+A)

  print(datetime.datetime.now(), ": Getting minor allele frequencies...")
  freqs = [min(x, 1-x) for x in freqs]

  # if only allele frequencies are desired, stop here
  print(datetime.datetime.now(), ": Saving results...")
  np.savetxt(output_prefix + "_freqs.csv", freqs, delimiter = ",", fmt='%s')

  # if its desired to estimate allele state by bayes theorem, go another step further

# functions to compare two hetmer tables

# Main function that collects everything together
@click.command(context_settings={'show_default': True})
#@click.option("-i", "--input", default=None, help="Path to input k-mer count table (KMC or Jellyfish)", multiple=False)
#@click.option("-t", "--type", default='kmer', type=click.Choice(['kmer', 'hetmer'], case_sensitive=False, help="Input type (het-mer count table or k-mer count table)", multiple=False)
@click.option('-k', '--kmer_table', help="Path to input k-mer count table (KMC or Jellyfish)")
@click.option('-t', '--hetmer_table', type=str, help="Prefix to hetmer tables (hetmers.py)")
@click.option('-y', '--analysis', type=click.Choice(['fst', 'freq'], case_sensitive=False), default = None, help="Type of analysis to do on hetmer table")
@click.option("-a", "--alleles", default=[2], help="Number of alleles", multiple=True)
@click.option("-m", "--minimum", default=5, help="Minimum required k-mer count", multiple=False)
@click.option("-o", "--output-prefix", help="Prefix for output files", required = True)
def main(kmer_table, hetmer_table, analysis, alleles, minimum, output_prefix):
  """
  Functions to find heterozygous pairs of k-mers (i.e. hetmers) in k-mer count tables and do some simple population genetics analyses.
  """
  # Make sure that only one analysis is being specified
  if (kmer_table is None and hetmer_table is None) or (kmer_table is not None and hetmer_table is not None):
    raise click.UsageError(f"Either 'kmer_table' or 'hetmer_table' must be specified, but not both.")

  if (hetmer_table is not None) and (analysis is None):
    raise click.UsageError(f"If 'hetmer_table' is specified, then 'analysis' must also be specified")

  if kmer_table is not None:
    get_hetmers(kmer_table, type, alleles, minimum, output_prefix)

  if hetmer_table is not None and (analysis == 'freq'):
    get_allele_freqs(hetmer_table, output_prefix)

  if hetmer_table is not None and (analysis == 'fst'):
    print("Working on it!")

if __name__ == '__main__':
    main()


