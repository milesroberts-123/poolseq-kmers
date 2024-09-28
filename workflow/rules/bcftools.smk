rule bcftools:
	input:
		"slim_{ID}.vcf"
	output:
		samplevcf = "samples_{ID}.vcf.gz",
		allelefreq = "slim_allele_freqs_{ID}.txt"
        threads: 1
        resources:
                mem_mb_per_cpu=8000,
                time=239
        conda:
                "../envs/bcftools.yml"
	shell:
		"""
		# compress simulation outputs
		bgzip slim_{wildcards.ID}.vcf
		tabix slim_{wildcards.ID}.vcf.gz

		# remove reference
		bcftools view --samples-file ../config/samples.txt -Oz -o {output.samplevcf} slim_{wildcards.ID}.vcf.gz

		# calculate allele frequencies
		bcftools +fill-tags {output.samplevcf} -Oz > samples_filled_{wildcards.ID}.vcf.gz
		bcftools query -f '%CHROM %POS %AF %AC\n' samples_filled_{wildcards.ID}.vcf.gz > {output.allelefreq}
		"""
