rule kmc_two_pop:
	input:
                pread1_p1 = "trimmed_paired_R1_{ID}_p1.fastq",
                pread2_p1 = "trimmed_paired_R2_{ID}_p1.fastq",
                uread1_p1 = "trimmed_unpaired_R1_{ID}_p1.fastq",
                uread2_p1 = "trimmed_unpaired_R2_{ID}_p1.fastq",
                pread1_p2 = "trimmed_paired_R1_{ID}_p2.fastq",
                pread2_p2 = "trimmed_paired_R2_{ID}_p2.fastq",
                uread1_p2 = "trimmed_unpaired_R1_{ID}_p2.fastq",
                uread2_p2 = "trimmed_unpaired_R2_{ID}_p2.fastq"
	output:
		p1="kmer_counts_{ID}_p1.txt",
		p2="kmer_counts_{ID}_p2.txt"
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
		kmc -ci1 -k31 {input.pread1_p1} tmp_R1_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}
		kmc -ci1 -k31 {input.pread2_p1} tmp_R2_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}
		kmc -ci1 -k31 {input.uread1_p1} tmp_u_R1_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}
		kmc -ci1 -k31 {input.uread2_p1} tmp_u_R2_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}

		kmc -ci1 -k31 {input.pread1_p2} tmp_R1_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}
		kmc -ci1 -k31 {input.pread2_p2} tmp_R2_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}
		kmc -ci1 -k31 {input.uread1_p2} tmp_u_R1_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}
		kmc -ci1 -k31 {input.uread2_p2} tmp_u_R2_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}

		# combine k-mer counts into one database
		kmc_tools simple tmp_R1_{wildcards.ID}_p1 tmp_R2_{wildcards.ID}_p1 union union_R1_R2_{wildcards.ID}_p1 &>> {log}
		kmc_tools simple union_R1_R2_{wildcards.ID}_p1 tmp_u_R1_{wildcards.ID}_p1 union union_R1_R2_u1_{wildcards.ID}_p1 &>> {log}
		kmc_tools simple union_R1_R2_u1_{wildcards.ID}_p1 tmp_u_R2_{wildcards.ID}_p1 union union_R1_R2_u1_u2_{wildcards.ID}_p1 &>> {log}

		kmc_tools simple tmp_R1_{wildcards.ID}_p2 tmp_R2_{wildcards.ID}_p2 union union_R1_R2_{wildcards.ID}_p2 &>> {log}
		kmc_tools simple union_R1_R2_{wildcards.ID}_p2 tmp_u_R1_{wildcards.ID}_p2 union union_R1_R2_u1_{wildcards.ID}_p2 &>> {log}
		kmc_tools simple union_R1_R2_u1_{wildcards.ID}_p2 tmp_u_R2_{wildcards.ID}_p2 union union_R1_R2_u1_u2_{wildcards.ID}_p2 &>> {log}

		# dump all k-mers to text file
		kmc_tools transform union_R1_R2_u1_u2_{wildcards.ID}_p1 dump {output.p1} &>> {log}
		
		kmc_tools transform union_R1_R2_u1_u2_{wildcards.ID}_p1 dump {output.p2} &>> {log}

		# delete tmp directory
		rm -r tmp_*_{wildcards.ID}_* union_*_{wildcards.ID}_*.kmc_*
		"""
