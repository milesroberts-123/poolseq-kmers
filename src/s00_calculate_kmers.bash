#
mamba activate bcftools
bgzip slim.vcf
tabix slim.vcf.gz

# remove reference
bcftools view --samples-file samples.txt -Oz -o samples.vcf.gz slim.vcf.gz

# calculate allele frequencies
bcftools +fill-tags samples.vcf.gz -Oz > samples_filled.vcf.gz
bcftools query -f '%CHROM %POS %AF %AC\n' samples_filled.vcf.gz > slim_allele_freqs.txt

## generate reads from slim simulation
#iss generate -g slim.fasta --cpus 4 --model miseq -n 300000 --abundance uniform --output reads
iss generate -g samples.fasta --cpus 8 --model miseq -n 166667 --abundance uniform --output reads

## trim reads
#fastp -u 40 -q 30 -l 31 -i reads_R1.fastq -I reads_R2.fastq -o trimmed_R1.fastq -O trimmed_R2.fastq
fastp -u 40 -q 30 -l 31 -i reads_R1.fastq -I reads_R2.fastq -o trimmed_paired_R1.fastq -O trimmed_paired_R2.fastq --unpaired-1 trimmed_unpaired_R1.fastq --unpaired-2 trimmed_unpaired_R2.fastq

## count k-mers
mkdir tmp_kmc

kmc -ci3 -k31 trimmed_paired_R1.fastq tmp_R1 tmp_kmc
kmc -ci3 -k31 trimmed_paired_R2.fastq tmp_R2 tmp_kmc
kmc -ci3 -k31 trimmed_unpaired_R1.fastq tmp_u_R1 tmp_kmc
kmc -ci3 -k31 trimmed_unpaired_R2.fastq tmp_u_R2 tmp_kmc

#kmc_tools simple tmp_R1 tmp_R2 union union_R1_R2
kmc_tools simple tmp_R1 tmp_R2 union union_R1_R2
kmc_tools simple union_R1_R2 tmp_u_R1 union union_R1_R2_u1
kmc_tools simple union_R1_R2_u1 tmp_u_R2 union union_R1_R2_u1_u2

kmc_tools transform union_R1_R2_u1_u2 dump kmer_counts.txt

# get counts of fake k-mers
#kmc -fm -ci1 -k31 slim.fasta tmp_fa tmp_kmc
#kmc_tools simple union_R1_R2 tmp_fa kmers_subtract fake_kmers
#kmc_tools transform fake_kmers dump fake_kmer_counts.txt

# convert tab to space delimiter
cat kmer_counts.txt | tr '\t' ' ' > kmer_counts.space

# loop over slim vcf
zcat samples_filled.vcf.gz | grep -v "^#" | cut -f 2 > snp_positions.txt

readarray posarray < snp_positions.txt

rm kmer_pairs.txt
for i in "${posarray[@]}"
do
echo Extracting k-mers for snp $i

center_start=$(($i-15))
center_end=$(($i+15))

cat samples.fasta | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u > center_kmer_pair.txt

join -t' ' -1 1 -2 1 -o 2.1,2.2 <(sort center_kmer_pair.txt) <(cat kmer_counts.space) > center_kmer_pair_counts.txt

echo $i $(cat center_kmer_pair_counts.txt | tr '\n' ' ') >> kmer_pairs.txt

done

# get key for k-mers to snp positions
readarray posarray < snp_positions.txt
rm center_kmer_pairs.txt
for i in "${posarray[@]}"
do
echo Extracting k-mers for snp $i

center_start=$(($i-15))
center_end=$(($i+15))

cat slim.fasta | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> center_kmer_pairs.txt
done

# run smudgeplot to extract pairs of k-mers that differ at the central bp
smudgeplot.py hetkmers --middle kmer_counts.txt