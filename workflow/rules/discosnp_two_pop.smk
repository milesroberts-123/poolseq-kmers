rule discosnp_two_pop:
	input:
                pread1_p1 = "trimmed_paired_R1_{ID}_p1.fastq",
                pread2_p1 = "trimmed_paired_R2_{ID}_p1.fastq",
                uread1_p1 = "trimmed_unpaired_R1_{ID}_p1.fastq",
                uread2_p1 = "trimmed_unpaired_R2_{ID}_p1.fastq",
                pread1_p2 = "trimmed_paired_R1_{ID}_p2.fastq",
                pread2_p2 = "trimmed_paired_R2_{ID}_p2.fastq",
                uread1_p2 = "trimmed_unpaired_R1_{ID}_p2.fastq",
                uread2_p2 = "trimmed_unpaired_R2_{ID}_p2.fastq",
		ref = "ref_{ID}_p1.fasta"
	output:
		tmpread_p1 = temp("tmp_read_set_{ID}_p1.fastq"),
		fasta_p1 = temp("discoRes_{ID}_p1_k_31_c_3_D_100_P_3_b_0_coherent.fa"),
		fof_p1 = temp("fof_{ID}_p1.txt"),
		vcf_p1 = "discoRes_{ID}_p1_k_31_c_3_D_100_P_3_b_0_coherent.vcf",
		tmpread_p2 = temp("tmp_read_set_{ID}_p2.fastq"),
		fasta_p2 = temp("discoRes_{ID}_p2_k_31_c_3_D_100_P_3_b_0_coherent.fa"),
		fof_p2 = temp("fof_{ID}_p2.txt"),
		vcf_p2 = "discoRes_{ID}_p2_k_31_c_3_D_100_P_3_b_0_coherent.vcf"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239,
		load = 1
	conda:
		"../envs/discosnp.yaml"
	log: 
		"logs/discosnp/{ID}.log"
	params:
		prefix_p1 = "discoRes_{ID}_p1",
		prefix_p2 = "discoRes_{ID}_p1"
	shell:
		"""
		# combine reads into one set
		cat {input.pread1_p1} {input.pread2_p1} {input.uread1_p1} {input.uread2_p1} > {output.tmpread_p1}
		cat {input.pread1_p2} {input.pread2_p2} {input.uread1_p2} {input.uread2_p2} > {output.tmpread_p2}

		# create file of files
		echo "{output.tmpread_p1}" > {output.fof_p1}
		echo "{output.tmpread_p2}" > {output.fof_p2}

		# run discosnp, with results for mapping SNPs to reference
		run_discoSnp++.sh -r {output.fof_p1} -c 3 -G {input.ref} -p {params.prefix_p1}
		run_discoSnp++.sh -r {output.fof_p2} -c 3 -G {input.ref} -p {params.prefix_p2}
		"""
