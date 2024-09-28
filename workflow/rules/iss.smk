rule iss:
	input:
		"samples_{ID}.fasta"
	output:
		"reads_{ID}_R1.fastq",
		"reads_{ID}_R2.fastq"
        threads: 1
        resources:
                mem_mb_per_cpu=8000,
                time=239
        conda:
                "../envs/iss.yml"
	shell:
		"iss generate -g {input} --cpus {threads} --model miseq -n 333334 --abundance uniform --output reads_{wildcards.ID}"
