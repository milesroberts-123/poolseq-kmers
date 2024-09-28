rule bcftools:
	input:
		"slim_{ID}.vcf",
	output:
		samplevcf = "samples_{ID}.vcf.gz",
		allelefreq = "slim_allele_freqs_{ID}.txt",
		filledvcf = "samples_filled_{ID}.vcf.gz"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/bcftools_{ID}.log"
	shell:
		"""
		# compress simulation outputs
		bgzip slim_{wildcards.ID}.vcf &> {log}
		tabix slim_{wildcards.ID}.vcf.gz &> {log}

		# remove reference
		bcftools view --samples-file ../config/ref.txt -Oz -o {output.samplevcf} slim_{wildcards.ID}.vcf.gz &> {log}

		# calculate allele frequencies
		bcftools +fill-tags {output.samplevcf} -Oz > {output.filledvcf} &> {log}
		bcftools query -f '%CHROM %POS %AF %AC\n' {output.filledvcf} > {output.allelefreq} &> {log}
		"""
