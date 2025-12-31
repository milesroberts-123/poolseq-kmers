#!/bin/bash

mamba activate bcftools

# path to reference genome
# chromosome names should be same as vcf
# all letters should be capital
ref=/global/scratch/users/milesroberts/tmp/arabidopsis_thaliana_tair10_ncbi/data/GCF_000001735.4/at_upper.fa

# loop through variant types from plantpan
#declare -a vartype=("snp" "ins" "cpg" "cpl" "del" "dup" "hdr" "inval" "invdp" "invtr" "tdm" "trans")
declare -a vartype=("snp" "ins" "del")

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

# realign indels and remove duplicatesi
echo Normalizing variants...
bcftools norm -m +any -d all -o norm.vcf.gz -O z -c s -f $ref all.vcf.gz

# filter out rare variants
echo Filtering out rare variants...
bcftools view -v snps,indels -o mac2.vcf -i 'MAC > 0' norm.vcf.gz

# remove any variants with complex mutations that overlap with the simpler mutations
#echo Removing overlapping variants...
#zcat mac2.vcf.gz | awk '$3 !~ /CPG|DUP|INV|TDM|CPL|TRANS/ || $1 ~ /^#/' > nonoverlap.vcf

# create genome file for bedtools
#echo Creating genome file...
#samtools faidx $ref
#awk '{print $1"\t"$2}' $ref.fai > $ref.genome

# create variant flanks
#echo Creating variant flanks...
#bcftools query -f '%CHROM\t%POS0\t%END\n' mac2.vcf.gz > mac2.bed
#bedtools flank -i mac2.bed -g at.genome -b 36 > flanks.bed

# subtract flanks from variants to get variants not near other variants
#echo Getting isolated variants...
#bedtools subtract -A -a mac2.vcf.gz -b flanks.bed > isolated.vcf
#bcftools view -h mac2.vcf.gz > header.vcf
#cat header.vcf isolated.vcf > isolated_with_header.vcf
#bgzip isolated_with_header.vcf
#tabix isolated_with_header.vcf.gz

# get only snps
#echo Extracting only snps...
#bcftools view -v snps isolated_with_header.vcf.gz -Oz -o snps_only.vcf.gz

# some non-snps were not filtered out still
#echo Removing DUP and TRANS...
#zgrep -v "DUP" snps_only.vcf.gz | grep -v "TRANS" > no_dup.vcf

# remove a snp near the chromosome end
# zgrep -v "SNP687663;SNP945145;SNP725163;SNP734774;SNP795013" no_dup.vcf.gz > no_ends.vcf

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

# realign indels and remove duplicates
#bcftools norm -d all -o final.vcf.gz -O z -c s -f at.fa capital.vcf.gz
