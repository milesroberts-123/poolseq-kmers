#!/bin/bash

mamba activate bcftools

# path to reference genome
# chromosome names should be same as vcf
# all letters should be capital
ref=/global/scratch/users/milesroberts/tmp/arabidopsis_thaliana_tair10_ncbi/data/GCF_000001735.4/at_upper.fa
#ref=/global/scratch/users/milesroberts/tmp/brassica_rapa_plantpan_data/br.genome.gene/Chiifu.genome.fasta
#ref=/global/scratch/users/milesroberts/tmp/brassica_napus_plantpan_data/bn.genome.gene/zs11_upper.fa
#ref=/global/scratch/users/milesroberts/tmp/solanum_lycopersicum_plantpan_data/to.genome.gene/SL5.genome.fasta

# loop through variant types from plantpan
#declare -a vartype=("snp" "ins" "cpg" "cpl" "del" "dup" "hdr" "inval" "invdp" "invtr" "tdm" "trans")
declare -a vartype=("snp" "ins" "del")
declare -a vartype=("cpg" "cpl")

for var in "${vartype[@]}"
do

echo $var

# Add GT fields
echo Adding GT field to vcf files
for FILE in *."$var".vcf
do
echo $FILE
SAMPLE=$(basename $FILE .$var.vcf)
echo $SAMPLE
sed "s:\tINFO:\tINFO\tFORMAT\t$SAMPLE:g" $FILE | sed "5i##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" > $(basename $FILE .vcf)_FORMAT.vcf
awk -F'\t' -v OFS="\t" '{ if(/^#/){print}else{$9="GT"$9;$10="1"$10;print}}' $(basename $FILE .vcf)_FORMAT.vcf > $(basename $FILE .vcf).GT.vcf
done

# bgzip files
echo bgziping vcf files...

for FILE in *"$var".GT.vcf
do
echo $FILE
bgzip $FILE
done

#index files
echo indexing files...

for FILE in *"$var".GT.vcf.gz
do
echo $FILE
tabix $FILE
done

# merge genotypes
echo Merging $var...
bcftools merge *$var.GT.vcf.gz --missing-to-ref -Oz -o all.$var.vcf.gz

done

for FILE in all.*.vcf.gz
do
echo $FILE
tabix $FILE
done

# concatenate all variants
echo Concatenating all variants...
bcftools concat -a -Oz -o all.vcf.gz all.*.vcf.gz

# realign indels and remove duplicates
# merge bi-alleleic records
echo Normalizing variants...
bcftools norm -m +any -d all -o norm.vcf.gz -c s -Oz -f $ref all.vcf.gz

# filter out rare variants
echo Filtering out rare variants...
bcftools view -v snps,indels -o mac2.vcf -i 'MAC > 0' norm.vcf.gz

# capitalize letters
echo Capitalize letters for alleles...
awk '{if ($0 ~ /^#/) print; else {split($4, a, ","); for (i=1; i<=length(a); i++) a[i]=toupper(a[i]); $4=a[1]; for (i=2; i<=length(a); i++) $4=$4","a[i]; split($5, b, ","); for (i=1; i<=length(b); i++) b[i]=toupper(b[i]); $5=b[1]; for (i=2; i<=length(b); i++) $5=$5","b[i]; print}}' mac2.vcf > final.vcf

# convert spaces to tabs
echo Convert spaces to tabs...
grep -v "^#" final.vcf | tr ' ' '\t' > body.vcf
grep "^#" final.vcf > header.vcf
cat header.vcf body.vcf > final.vcf
bgzip final.vcf
tabix final.vcf.gz

