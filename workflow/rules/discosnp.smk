rule discosnp:
	input:
                pread1 = "trimmed_paired_R1_{ID}.fastq",
                pread2 = "trimmed_paired_R2_{ID}.fastq",
                uread1 = "trimmed_unpaired_R1_{ID}.fastq",
                uread2 = "trimmed_unpaired_R2_{ID}.fastq",
		ref = "ref_{ID}.fasta"
	output:
		tmpread = temp("tmp_read_set_{ID}.fastq"),
		fasta = temp("discoRes_{ID}_k_31_c_3_D_100_P_3_b_0_coherent.fa"),
		fof = temp("fof_{ID}.txt"),
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
	params:
		prefix = "discoRes_{ID}"
	shell:
		"""
		# combine reads into one set
		cat {input.pread1} {input.pread2} {input.uread1} {input.uread2} > {output.tmpread}

		# create file of files
		echo "{output.tmpread}" > {output.fof}

		# run discosnp, with results for mapping SNPs to reference
		run_discoSnp++.sh -r {output.fof} -c 3 -G {input.ref} -p {params.prefix}
		"""
