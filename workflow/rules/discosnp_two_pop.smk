rule discosnp_two_pop:
	input:
                pread1_p1 = "fastp_results/trimmed_paired_R1_{ID}_p1.fastq",
                pread2_p1 = "fastp_results/trimmed_paired_R2_{ID}_p1.fastq",
                uread1_p1 = "fastp_results/trimmed_unpaired_R1_{ID}_p1.fastq",
                uread2_p1 = "fastp_results/trimmed_unpaired_R2_{ID}_p1.fastq",
                pread1_p2 = "fastp_results/trimmed_paired_R1_{ID}_p2.fastq",
                pread2_p2 = "fastp_results/trimmed_paired_R2_{ID}_p2.fastq",
                uread1_p2 = "fastp_results/trimmed_unpaired_R1_{ID}_p2.fastq",
                uread2_p2 = "fastp_results/trimmed_unpaired_R2_{ID}_p2.fastq",
		ref = "seqkit_results/ref_{ID}_p1.fasta",
		amb = "seqkit_results/ref_{ID}_p1.fasta.amb",
		ann = "seqkit_results/ref_{ID}_p1.fasta.ann",
		bwt = "seqkit_results/ref_{ID}_p1.fasta.bwt",
		pac = "seqkit_results/ref_{ID}_p1.fasta.pac",
		sa = "seqkit_results/ref_{ID}_p1.fasta.sa"
	output:
		#tmpread_p1 = "tmp_read_set_{ID}_p1.fastq",
		fasta_p1 = temp("discoRes_{ID}_p1_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.fa"),
		fof = temp("fof_{ID}.txt"),
		fof_p1 = temp("fof_{ID}_p1.txt"),
		vcf_p1 = "discoRes_{ID}_p1_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.vcf",
		#tmpread_p2 = "tmp_read_set_{ID}_p2.fastq",
		fasta_p2 = temp("discoRes_{ID}_p2_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.fa"),
		fof_p2 = temp("fof_{ID}_p2.txt"),
		vcf_p2 = "discoRes_{ID}_p2_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.vcf"
	threads: 1
	resources:
		mem_mb_per_cpu=16000,
		time=239,
		load = 1
	conda:
		"../envs/discosnp.yaml"
	log: 
		"logs/discosnp/{ID}.log"
	params:
		prefix_p1 = "discoRes_{ID}_p1",
		prefix_p2 = "discoRes_{ID}_p2",
		mincount = config["mincount"],
		k = config["k"]
	priority: 100
	shell:
		"""
		# create file of files
		echo "{output.fof_p1}" > {output.fof} 

		# check each file for being empty, use only non-empty files
		if [ -s {input.pread1_p1} ]; then
			echo {input.pread1_p1} >> {output.fof_p1}
		fi

		if [ -s {input.pread2_p1} ]; then
			echo {input.pread2_p1} >> {output.fof_p1}
		fi

		if [ -s {input.uread1_p1} ]; then
			echo {input.uread1_p1} >> {output.fof_p1}
		fi

		if [ -s {input.uread2_p1} ]; then
			echo {input.uread2_p1} >> {output.fof_p1}
		fi

		# run discosnp, with results for mapping SNPs to reference
		run_discoSnp++.sh -r {output.fof} -c {params.mincount} -G {input.ref} -p {params.prefix_p1}

		# repeat for population 2
		echo "{output.fof_p2}" > {output.fof}

		# check each file for being empty, use only non-empty files
		if [ -s {input.pread1_p2} ]; then
			echo {input.pread1_p2} >> {output.fof_p2}
		fi

		if [ -s {input.pread2_p2} ]; then
			echo {input.pread2_p2} >> {output.fof_p2}
		fi

		if [ -s {input.uread1_p2} ]; then
			echo {input.uread1_p2} >> {output.fof_p2}
		fi

		if [ -s {input.uread2_p2} ]; then
			echo {input.uread2_p2} >> {output.fof_p2}
		fi

		run_discoSnp++.sh -r {output.fof} -c {params.mincount} -k {params.k} -G {input.ref} -p {params.prefix_p2}
		"""
