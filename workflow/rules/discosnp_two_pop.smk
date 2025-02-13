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
		ref = "ref_{ID}_p1.fasta",
		amb = "ref_{ID}_p1.fasta.amb",
		ann = "ref_{ID}_p1.fasta.ann",
		bwt = "ref_{ID}_p1.fasta.bwt",
		pac = "ref_{ID}_p1.fasta.pac",
		sa = "ref_{ID}_p1.fasta.sa"
	output:
		#tmpread_p1 = "tmp_read_set_{ID}_p1.fastq",
		fasta_p1 = temp("discoRes_{ID}_p1_k_31_c_3_D_100_P_3_b_0_coherent.fa"),
		fof = "fof_{ID}.txt",
		fof_p1 = "fof_{ID}_p1.txt",
		vcf_p1 = "discoRes_{ID}_p1_k_31_c_3_D_100_P_3_b_0_coherent.vcf",
		#tmpread_p2 = "tmp_read_set_{ID}_p2.fastq",
		fasta_p2 = temp("discoRes_{ID}_p2_k_31_c_3_D_100_P_3_b_0_coherent.fa"),
		fof_p2 = "fof_{ID}_p2.txt",
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
		prefix_p2 = "discoRes_{ID}_p2"
	shell:
		"""
		# create file of files
		echo "{output.fof_p1}" > {output.fof} 
		echo -e "{input.pread1_p1}\n{input.pread2_p1}\n{input.uread1_p1}\n{input.uread2_p1}" > {output.fof_p1}

		# run discosnp, with results for mapping SNPs to reference
		run_discoSnp++.sh -r {output.fof} -c 3 -G {input.ref} -p {params.prefix_p1}

		# repeat for population 2
		echo "{output.fof_p2}" > {output.fof}
		echo -e "{input.pread1_p2}\n{input.pread2_p2}\n{input.uread1_p2}\n{input.uread2_p2}" > {output.fof_p2}
		run_discoSnp++.sh -r {output.fof} -c 3 -G {input.ref} -p {params.prefix_p2}
		"""
