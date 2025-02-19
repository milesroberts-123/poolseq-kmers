rule discosnp_one_pop:
	input:
                pread1 = "trimmed_paired_R1_{ID}.fastq",
                pread2 = "trimmed_paired_R2_{ID}.fastq",
                uread1 = "trimmed_unpaired_R1_{ID}.fastq",
                uread2 = "trimmed_unpaired_R2_{ID}.fastq",
		ref = "ref_{ID}.fasta",
		amb = "ref_{ID}.fasta.amb",
		ann = "ref_{ID}.fasta.ann",
		bwt = "ref_{ID}.fasta.bwt",
		pac = "ref_{ID}.fasta.pac",
		sa = "ref_{ID}.fasta.sa"
	output:
		#tmpread = temp("tmp_read_set_{ID}.fastq"),
		fasta = temp("discoRes_{ID}_k_31_c_3_D_100_P_3_b_0_coherent.fa"),
		fof = temp("fof_{ID}.txt"),
		fof_reads = temp("fof_reads_{ID}.txt"),
		vcf = "discoRes_{ID}_k_31_c_3_D_100_P_3_b_0_coherent.vcf"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239,
		load = 1
	conda:
		"../envs/discosnp.yaml"
	log: 
		"logs/discosnp/{ID}.log"
	benchmark:
		"benchmarks/discosnp/{ID}.bench"
	params:
		prefix = "discoRes_{ID}"
	shell:
		"""
		# create file of files
		echo "{output.fof_reads}" > {output.fof}
		echo -e "{input.pread1}\n{input.pread2}\n{input.uread1}\n{input.uread2}" > {output.fof_reads}

		# run discosnp, with results for mapping SNPs to reference
		run_discoSnp++.sh -r {output.fof} -c 3 -G {input.ref} -p {params.prefix}
		"""
