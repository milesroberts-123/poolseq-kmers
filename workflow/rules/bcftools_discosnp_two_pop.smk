rule bcftools_discosnp_two_pop:
	input:
		ref = "ref_{ID}.fasta",
		vcf_p1 = "discoRes_{ID}_p1_k_31_c_3_D_100_P_3_b_0_coherent.vcf",
		vcf_p2 = "discoRes_{ID}_p2_k_31_c_3_D_100_P_3_b_0_coherent.vcf",
	output:
		fai = temp("ref_{ID}.fasta.fai"),
		header_p1 = temp("discoRes_header_{ID}_p1.vcf"),
		header_p2 = temp("discoRes_header_{ID}_p2.vcf"),
		bgzip_p1 = temp("discoRes_sorted_{ID}_p1.vcf.gz"),
		bgzip_p2 = temp("discoRes_sorted_{ID}_p2.vcf.gz"),
		tbi_p1 = temp("discoRes_sorted_{ID}_p1.vcf.gz.tbi"),
		tbi_p2 = temp("discoRes_sorted_{ID}_p2.vcf.gz.tbi"),
		final_p1 = "discoRes_ad_{ID}_p1.txt",
		final_p2 = "discoRes_ad_{ID}_p2.txt"
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
		bcftools reheader --fai {output.fai} {input.vcf_p1} | sed 's:TySNP:Ty=SNP:g' > {output.header_p1}
		bcftools reheader --fai {output.fai} {input.vcf_p2} | sed 's:TySNP:Ty=SNP:g' > {output.header_p2}

		# sort vcf
		bcftools sort {output.header_p1} > discoRes_sorted_{wildcards.ID}_p1.vcf
		bcftools sort {output.header_p2} > discoRes_sorted_{wildcards.ID}_p2.vcf

		# index vcf
		bgzip discoRes_sorted_{wildcards.ID}_p1.vcf
		tabix {output.bgzip_p1}

		bgzip discoRes_sorted_{wildcards.ID}_p2.vcf
		tabix {output.bgzip_p2}

		# output allele depths
		bcftools query -f '%CHROM %POS [ %AD]\n' {output.bgzip_p1} | sed 's:,:\t:g' > {output.final_p1}
		bcftools query -f '%CHROM %POS [ %AD]\n' {output.bgzip_p2} | sed 's:,:\t:g' > {output.final_p2}
		"""
