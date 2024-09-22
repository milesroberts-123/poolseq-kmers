#
mamba activate bcftools
bgzip slim.vcf
tabix slim.vcf
bcftools +fill-tags slim.vcf.gz > slim_filled.vcf.gz
bcftools query -f '%CHROM %POS %AF\n' slim_filled.vcf.gz > slim_allele_freqs.txt

## generate reads from slim simulation
iss generate -g slim.fasta --cpus 4 --model miseq -n 300000 --abundance uniform --output reads

## trim reads
fastp -u 40 -q 30 -l 31 -i reads_R1.fastq -I reads_R2.fastq -o trimmed_R1.fastq -O trimmed_R2.fastq

## count k-mers
mkdir tmp_kmc

kmc -ci3 -k31 trimmed_R1.fastq tmp_R1 tmp_kmc
kmc -ci3 -k31 trimmed_R2.fastq tmp_R2 tmp_kmc
kmc_tools simple tmp_R1 tmp_R2 union union_R1_R2
kmc_tools transform union_R1_R2 dump kmer_counts.txt

# get counts of fake k-mers
#kmc -fm -ci1 -k31 slim.fasta tmp_fa tmp_kmc
#kmc_tools simple union_R1_R2 tmp_fa kmers_subtract fake_kmers
#kmc_tools transform fake_kmers dump fake_kmer_counts.txt

# convert tab to space delimiter
cat kmer_counts.txt | tr '\t' ' ' > kmer_counts.space

# loop over slim vcf
zcat slim.vcf.gz | grep -v "^#" | cut -f 2 > snp_positions.txt

readarray posarray < snp_positions.txt

rm kmer_pairs.txt
for i in "${posarray[@]}"
do
echo Extracting k-mers for snp $i

center_start=$(($i-15))
center_end=$(($i+15))

cat slim.fasta | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u > center_kmer_pair.txt

join -t' ' -1 1 -2 1 -o 2.1,2.2 <(sort center_kmer_pair.txt) <(cat kmer_counts.space) > center_kmer_pair_counts.txt

echo $i $(cat center_kmer_pair_counts.txt | tr '\n' ' ') >> kmer_pairs.txt

done
