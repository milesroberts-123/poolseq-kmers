rule kmc:
	input:
                pread1 = "trimmed_paired_R1_{ID}.fastq",
                pread2 = "trimmed_paired_R2_{ID}.fastq",
                uread1 = "trimmed_unpaired_R1_{ID}.fastq",
                uread2 = "trimmed_unpaired_R2_{ID}.fastq"
	output:
		"kmer_counts_{ID}.txt"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/kmc.yaml"
	log: 
		"logs/kmc/{ID}.log"
	shell:
		"""
		# create directory
		if [ -d "tmp_kmc_{wildcards.ID}" ]; then
			rm -r tmp_kmc_{wildcards.ID}
		fi

		mkdir tmp_kmc_{wildcards.ID}

		# count k-mers
		kmc -ci1 -k31 {input.pread1} tmp_R1_{wildcards.ID} tmp_kmc_{wildcards.ID} &>> {log}
		kmc -ci1 -k31 {input.pread2} tmp_R2_{wildcards.ID} tmp_kmc_{wildcards.ID} &>> {log}
		kmc -ci1 -k31 {input.uread1} tmp_u_R1_{wildcards.ID} tmp_kmc_{wildcards.ID} &>> {log}
		kmc -ci1 -k31 {input.uread2} tmp_u_R2_{wildcards.ID} tmp_kmc_{wildcards.ID} &>> {log}

		# combine k-mer counts into one database
		kmc_tools simple tmp_R1_{wildcards.ID} tmp_R2_{wildcards.ID} union union_R1_R2_{wildcards.ID} &>> {log}
		kmc_tools simple union_R1_R2_{wildcards.ID} tmp_u_R1_{wildcards.ID} union union_R1_R2_u1_{wildcards.ID} &>> {log}
		kmc_tools simple union_R1_R2_u1_{wildcards.ID} tmp_u_R2_{wildcards.ID} union union_R1_R2_u1_u2_{wildcards.ID} &>> {log}

		# dump all k-mers to text file
		kmc_tools transform union_R1_R2_u1_u2_{wildcards.ID} dump {output} &>> {log}

		# delete tmp directory
		rm -r tmp_*_{wildcards.ID}
		"""
