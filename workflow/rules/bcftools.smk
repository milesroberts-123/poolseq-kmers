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
	log:
		"logs/bcftools_{ID}.log"
	shell:
		"""
		# compress simulation outputs
		bgzip slim_{wildcards.ID}.vcf &> {log}
		tabix slim_{wildcards.ID}.vcf.gz &> {log}

		# remove reference
		bcftools view --samples-file ../config/samples.txt -Oz -o {output.samplevcf} slim_{wildcards.ID}.vcf.gz &> {log}

		# calculate allele frequencies
		bcftools +fill-tags {output.samplevcf} -Oz > samples_filled_{wildcards.ID}.vcf.gz &> {log}
		bcftools query -f '%CHROM %POS %AF %AC\n' samples_filled_{wildcards.ID}.vcf.gz > {output.allelefreq} &> {log}
		"""
