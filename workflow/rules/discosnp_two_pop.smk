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
		tmpread_p1 = temp("discoRes_{ID}_p1_read_files_correspondance.txt"),
		fasta_p1 = temp("discoRes_{ID}_p1_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.fa"),
		un_fasta_p1 = temp("discoRes_{ID}_p1_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_uncoherent.fa"),
		fof1 = temp("fof1_{ID}.txt"),
		fof2 = temp("fof2_{ID}.txt"),
		fof_p1 = temp("fof_{ID}_p1.txt"),
		vcf_p1 = temp("discoRes_{ID}_p1_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.vcf"),
		igv_vcf_p1 = temp("discoRes_{ID}_p1_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent_for_IGV.vcf"),
		sam_p1 = temp("discoRes_{ID}_p1_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherentBWA_MEM.sam"),
		cov_h5_p1 = temp("discoRes_{ID}_p1_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_cov.h5"),
		tmpread_p2= temp("discoRes_{ID}_p2read_files_correspondance.txt"),
		fasta_p2 = temp("discoRes_{ID}_p2_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.fa"),
		un_fasta_p2 = temp("discoRes_{ID}_p2_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_uncoherent.fa"),
		fof_p2 = temp("fof_{ID}_p2.txt"),
		vcf_p2 = temp("discoRes_{ID}_p2_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.vcf"),
		igv_vcf_p2 = temp("discoRes_{ID}_p2_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent_for_IGV.vcf"),
		sam_p2 = temp("discoRes_{ID}_p2_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherentBWA_MEM.sam"),
		cov_h5_p2 = temp("discoRes_{ID}_p2_k_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_cov.h5")
	threads: 1
	resources:
		mem_mb_per_cpu=16000,
		time=239,
		#load = 1
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
		# create temp directory
		if [ -d "tmp_discosnp_{wildcards.ID}" ]; then
			rm -r tmp_discosnp_{wildcards.ID}
		fi

		mkdir tmp_discosnp_{wildcards.ID}

		# POPULATION ONE

		# create file of files
		echo "../{output.fof_p1}" > {output.fof1} 

		# check each file for being empty, use only non-empty files
		if [ -s {input.pread1_p1} ]; then
			echo "../{input.pread1_p1}" >> {output.fof_p1}
		fi

		if [ -s {input.pread2_p1} ]; then
			echo "../{input.pread2_p1}" >> {output.fof_p1}
		fi

		if [ -s {input.uread1_p1} ]; then
			echo "../{input.uread1_p1}" >> {output.fof_p1}
		fi

		if [ -s {input.uread2_p1} ]; then
			echo "../{input.uread2_p1}" >> {output.fof_p1}
		fi

		# POPULATION 2

		# create file of files
		echo "../{output.fof_p2}" > {output.fof2}

		# check each file for being empty, use only non-empty files
		if [ -s {input.pread1_p2} ]; then
			echo "../{input.pread1_p2}" >> {output.fof_p2}
		fi

		if [ -s {input.pread2_p2} ]; then
			echo "../{input.pread2_p2}" >> {output.fof_p2}
		fi

		if [ -s {input.uread1_p2} ]; then
			echo "../{input.uread1_p2}" >> {output.fof_p2}
		fi

		if [ -s {input.uread2_p2} ]; then
			echo "../{input.uread2_p2}" >> {output.fof_p2}
		fi


		# run discosnp, with results for mapping SNPs to reference
		cd tmp_discosnp_{wildcards.ID}

		run_discoSnp++.sh -r ../{output.fof1} -c {params.mincount} -k {params.k} -G ../{input.ref} -p {params.prefix_p1} &>> ../{log}

		run_discoSnp++.sh -r ../{output.fof2} -c {params.mincount} -k {params.k} -G ../{input.ref} -p {params.prefix_p2} &>> ../{log}

		# move output from temp directory
		mv {params.prefix_p1}* ..
		mv {params.prefix_p2}* ..

		# delete temporary directory
		cd ..
		rm -r tmp_discosnp_{wildcards.ID}
		"""
