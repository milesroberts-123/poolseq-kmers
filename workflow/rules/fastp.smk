rule fastp:
	input:
		read1 = "reads_{ID}_R1.fastq",
		read2 = "reads_{ID}_R2.fastq"
	output:
		pread1 = "trimmed_paired_R1_{ID}.fastq",
		pread2 = "trimmed_paired_R2_{ID}.fastq",
		uread1 = "trimmed_unpaired_R1_{ID}.fastq",
		uread2 = "trimmed_unpaired_R2_{ID}.fastq"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/fastp.yaml"
	log: 
		"logs/fastp_{ID}.log"
	shell:
		"fastp -u 40 -q 30 -l 31 -i {input.read1} -I {input.read2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2} &> {log}"
