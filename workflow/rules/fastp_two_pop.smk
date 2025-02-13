rule fastp_two_pop:
	input:
		read1_p1 = "reads_{ID}_p1_R1.fastq",
		read2_p1 = "reads_{ID}_p1_R2.fastq",
		read1_p2 = "reads_{ID}_p2_R1.fastq",
		read2_p2 = "reads_{ID}_p2_R2.fastq"
	output:
		pread1_p1 = temp("trimmed_paired_R1_{ID}_p1.fastq"),
		pread2_p1 = temp("trimmed_paired_R2_{ID}_p1.fastq"),
		uread1_p1 = temp("trimmed_unpaired_R1_{ID}_p1.fastq"),
		uread2_p1 = temp("trimmed_unpaired_R2_{ID}_p1.fastq"),
		pread1_p2 = temp("trimmed_paired_R1_{ID}_p2.fastq"),
		pread2_p2 = temp("trimmed_paired_R2_{ID}_p2.fastq"),
		uread1_p2 = temp("trimmed_unpaired_R1_{ID}_p2.fastq"),
		uread2_p2 = temp("trimmed_unpaired_R2_{ID}_p2.fastq"),
		json_p1 = "{ID}_p1.json",
		json_p2 = "{ID}_p2.json"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/fastp.yaml"
	log: 
		"logs/fastp/{ID}.log"
	shell:
		"""
		fastp -u 40 -q 30 -l 31 --dedup --correction --json {output.json_p1} -i {input.read1_p1} -I {input.read2_p1} -o {output.pread1_p1} -O {output.pread2_p1} --unpaired1 {output.uread1_p1} --unpaired2 {output.uread2_p1} &> {log}

		fastp -u 40 -q 30 -l 31 --dedup --correction --json {output.json_p2} -i {input.read1_p2} -I {input.read2_p2} -o {output.pread1_p2} -O {output.pread2_p2} --unpaired1 {output.uread1_p2} --unpaired2 {output.uread2_p2} &> {log}
		"""
