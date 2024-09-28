rule compress:
	input:
		"slim_{ID}.vcf",
	output:
		"slim_{ID}.vcf.gz",
		"slim_{ID}.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/compress/{ID}.log"
	shell:
		"""
			bgzip {input} &> {log}
			tabix {input}.gz &> {log}
		"""
