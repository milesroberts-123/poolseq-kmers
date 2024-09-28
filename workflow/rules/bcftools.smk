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
		"logs/bcftools/{ID}.log"
	shell:
		"""
		# compress simulation outputs
		if [ -e "{input}.gz" ]; then
			bgzip {input} &> {log}
		fi

		if [ -e "{input}.gz.tbi" ]; then
			tabix {input}.gz &> {log}
		fi

		# remove reference
		bcftools view --samples-file ../config/ref.txt -Oz -o {output.samplevcf} {input}.gz &> {log}

		# calculate allele frequencies
		bcftools +fill-tags {output.samplevcf} -Oz > {output.filledvcf} &> {log}
		bcftools query -f '%CHROM %POS %AF %AC\n' {output.filledvcf} > {output.allelefreq} &> {log}
		"""
