rule bcftools_discosnp:
	input:
		ref = "ref_{ID}.fasta",
		vcf = "discoRes_{ID}_k_31_c_3_D_100_P_3_b_0_coherent.vcf",
	output:
		fai = temp("ref_{ID}.fasta.fai"),
		header = temp("discoRes_header_{ID}.vcf"),
		sorted = temp("discoRes_sorted_{ID}.vcf"),
		bgzip = temp("discoRes_sorted_{ID}.vcf.gz"),
		tbi = temp("discoRes_sorted_{ID}.vcf.gz.tbi"),
		final = "discoRes_ad_{ID}.txt"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/bcftools_discosnp/{ID}.log"
	shell:
		"""
		# index reference
		samtools faidx {input.ref}

		# reheader
		# also, fix type definition for SNP < k bp apart
		bcftools reheader --fai {output.fai} {input.vcf} | sed 's:TySNP:Ty=SNP:g' > {output.header}

		# sort vcf
		bcftools sort {output.header} > {output.sorted}

		# index vcf
		bgzip {output.sorted}
		tabix {output.bgzip}

		# output allele depths
		bcftools query -f '%CHROM %POS [ %AD]\n' {output.bgzip} | sed 's:,:\t:g' > {output.final}
		"""
