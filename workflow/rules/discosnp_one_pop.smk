rule discosnp_one_pop:
	input:
                pread1 = "fastp_results/trimmed_paired_R1_{ID}.fastq",
                pread2 = "fastp_results/trimmed_paired_R2_{ID}.fastq",
                uread1 = "fastp_results/trimmed_unpaired_R1_{ID}.fastq",
                uread2 = "fastp_results/trimmed_unpaired_R2_{ID}.fastq",
		ref = "seqkit_results/ref_{ID}.fasta",
		amb = "seqkit_results/ref_{ID}.fasta.amb",
		ann = "seqkit_results/ref_{ID}.fasta.ann",
		bwt = "seqkit_results/ref_{ID}.fasta.bwt",
		pac = "seqkit_results/ref_{ID}.fasta.pac",
		sa = "seqkit_results/ref_{ID}.fasta.sa"
	output:
		#tmpread = temp("tmp_read_set_{ID}.fastq"),
		fasta = temp("discoRes_{ID}_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.fa"),
		fof = temp("fof_{ID}.txt"),
		fof_reads = temp("fof_reads_{ID}.txt"),
		vcf = "discoRes_{ID}_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent.vcf",
		h5 = temp("discoRes_{ID}_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_cov.h5"),
		sam = temp("discoRes_{ID}_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherentBWA_MEM.sam"),
		igv = temp("discoRes_{ID}_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_coherent_for_IGV.vcf"),
		uncofa = temp("discoRes_{ID}_" + str(config["k"]) + "_c_" + str(config["mincount"]) + "_D_100_P_3_b_0_uncoherent.fa")
	threads: 1
	resources:
		mem_mb_per_cpu=16000,
		time=239,
		load = 1
	conda:
		"../envs/discosnp.yaml"
	log: 
		"logs/discosnp/{ID}.log"
	benchmark:
		"benchmarks/discosnp/{ID}.bench"
	params:
		prefix = "discoRes_{ID}",
		mincount = config["mincount"],
		k = config["k"]
	priority: 100
	shell:
		"""
		# create file of files
		echo "{output.fof_reads}" > {output.fof}

		# check each file for being empty, use only non-empty files
		if [ -s {input.pread1} ]; then
			echo {input.pread1} >> {output.fof_reads}
		fi

		if [ -s {input.pread2} ]; then
			echo {input.pread2} >> {output.fof_reads}
		fi

		if [ -s {input.uread1} ]; then
			echo {input.uread1} >> {output.fof_reads}
		fi

		if [ -s {input.uread2} ]; then
			echo {input.uread2} >> {output.fof_reads}
		fi

		# run discosnp, with results for mapping SNPs to reference
		run_discoSnp++.sh -r {output.fof} -c {params.mincount} -k {params.k} -G {input.ref} -p {params.prefix} &> {log}
		"""
