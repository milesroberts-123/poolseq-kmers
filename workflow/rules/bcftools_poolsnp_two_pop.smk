rule bcftools_poolsnp_two_pop:
	input:
		vcf_p1 = "{ID}_p1_poolsnp_output.vcf.gz",
		vcf_p2 = "{ID}_p2_poolsnp_output.vcf.gz"
	output:
		tbi_p1 = temp("{ID}_p1_poolsnp_output.vcf.gz.tbi"),
		tbi_p2 = temp("{ID}_p2_poolsnp_output.vcf.gz.tbi"),
		final_p1 = "poolsnp_final_{ID}_p1.txt",
		final_p2 = "poolsnp_final_{ID}_p2.txt"
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
		tabix {input.vcf_p1}
		tabix {input.vcf_p2}

		# output allele depths
		bcftools query -f '%CHROM %POS [ %AD]\n' {input.vcf_p1} | sed 's:,:\t:g' > {output.final_p1}
		bcftools query -f '%CHROM %POS [ %AD]\n' {input.vcf_p2} | sed 's:,:\t:g' > {output.final_p2}
		"""
