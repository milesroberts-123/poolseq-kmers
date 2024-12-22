def get_samples(wildcards):
	# get sample size
        n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
	n = int(n.iloc[0])

	# create list of sample names from slim convention
	samples = list(range(1, n + 1))
	samples = ["i" + str(x) for x in samples] 

	# create comma-sep list for bcftools
	samples = ','.join(samples)

        return samples

rule bcftools_slim_two_pop:
	input:
		compvcf = "slim_{ID}.vcf.gz",
		vcfidx = "slim_{ID}.vcf.gz.tbi",
	output:
		allelefreq_p1 = "slim_allele_freqs_{ID}_p1.txt",
		allelefreq_p2 = "slim_allele_freqs_{ID}_p2.txt",
		filledvcf = temp("samples_filled_{ID}_p1p2.vcf.gz")
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/bcftools_slim/{ID}.log"
	params:
		samples=get_samples
	shell:
		"""
		# calculate allele frequencies
		bcftools +fill-tags {input.compvcf} -Oz -o {output.filledvcf} &> {log}

		# split allele frequencies by populations
		bcftools query -s {params.samples} -f '%CHROM %POS %NS %AF %AC\n' -o {output.allelefreq_p1} {output.filledvcf} &> {log}
		bcftools query -s ^{params.samples} -f '%CHROM %POS %NS %AF %AC\n' -o {output.allelefreq_p2} {output.filledvcf} &> {log}
		"""
