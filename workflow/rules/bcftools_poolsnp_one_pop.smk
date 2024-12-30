rule bcftools_poolsnp_one_pop:
	input:
		ref = "ref_{ID}.fasta",
		vcf = "{ID}_poolsnp_output.vcf.gz",
	output:
		tbi = temp("{ID}_poolsnp_output.vcf.gz.tbi"),
		final = "poolsnp_final_{ID}.txt"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/bcftools_poolsnp/{ID}.log"
	shell:
		"""
		# index vcf
		tabix {input.vcf}

		# output allele depths
		bcftools query -f '%CHROM %POS [ %AD]\n' {output.bgzip} | sed 's:,:\t:g' > {output.final}
		"""
