rule bcftools_slim:
	input:
		compvcf = "slim_{ID}.vcf.gz",
		vcfidx = "slim_{ID}.vcf.gz.tbi",
	output:
		samplevcf = temp("samples_{ID}.vcf.gz"),
		allelefreq = "slim_allele_freqs_{ID}.txt",
		filledvcf = temp("samples_filled_{ID}.vcf.gz")
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/bcftools_slim/{ID}.log"
	shell:
		"""
		# remove reference
		bcftools view --samples-file ^../config/ref.txt -Oz -o {output.samplevcf} {input.compvcf} &> {log}

		# calculate allele frequencies
		bcftools +fill-tags {output.samplevcf} -Oz -o {output.filledvcf} &> {log}
		bcftools query -f '%CHROM %POS %AF %AC\n' -o {output.allelefreq} {output.filledvcf} &> {log}
		"""
