rule fastp_two_pop:
	input:
		read1_p1 = "iss_results/reads_{ID}_p1_R1.fastq",
		read2_p1 = "iss_results/reads_{ID}_p1_R2.fastq",
		read1_p2 = "iss_results/reads_{ID}_p2_R1.fastq",
		read2_p2 = "iss_results/reads_{ID}_p2_R2.fastq"
	output:
		dpread1_p1 = temp("fastp_results/dedup_paired_R1_{ID}_p1.fastq"),
		dpread2_p1 = temp("fastp_results/dedup_paired_R2_{ID}_p1.fastq"),
		duread1_p1 = temp("fastp_results/dedup_unpaired_R1_{ID}_p1.fastq"),
		duread2_p1 = temp("fastp_results/dedup_unpaired_R2_{ID}_p1.fastq"),
		dpread1_p2 = temp("fastp_results/dedup_paired_R1_{ID}_p2.fastq"),
		dpread2_p2 = temp("fastp_results/dedup_paired_R2_{ID}_p2.fastq"),
		duread1_p2 = temp("fastp_results/dedup_unpaired_R1_{ID}_p2.fastq"),
		duread2_p2 = temp("fastp_results/dedup_unpaired_R2_{ID}_p2.fastq"),
		pread1_p1 = temp("fastp_results/trimmed_paired_R1_{ID}_p1.fastq"),
		pread2_p1 = temp("fastp_results/trimmed_paired_R2_{ID}_p1.fastq"),
		uread1_p1 = temp("fastp_results/trimmed_unpaired_R1_{ID}_p1.fastq"),
		uread2_p1 = temp("fastp_results/trimmed_unpaired_R2_{ID}_p1.fastq"),
		pread1_p2 = temp("fastp_results/trimmed_paired_R1_{ID}_p2.fastq"),
		pread2_p2 = temp("fastp_results/trimmed_paired_R2_{ID}_p2.fastq"),
		uread1_p2 = temp("fastp_results/trimmed_unpaired_R1_{ID}_p2.fastq"),
		uread2_p2 = temp("fastp_results/trimmed_unpaired_R2_{ID}_p2.fastq"),
		jsonR1R2_p1 = "fastp_results/{ID}_R1R2_p1.json",
		jsonU1_p1 = "fastp_results/{ID}_U1_p1.json",
		jsonU2_p1 = "fastp_results/{ID}_U2_p1.json",
		jsonR1R2_p2 = "fastp_results/{ID}_R1R2_p2.json",
		jsonU1_p2 = "fastp_results/{ID}_U1_p2.json",
		jsonU2_p2 = "fastp_results/{ID}_U2_p2.json",
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/fastp.yaml"
	log: 
		"logs/fastp/{ID}.log"
	params:
		unqualLimit = config["unqualLimit"],
		k = config["k"],
		qualThresh = config["qualThresh"],
		windowLength = config["windowLength"],
	shell:
		"""
		# population one

		## deduplicate and correct
		fastp -u {params.unqualLimit} -q {params.qualThresh} --dedup --correction -i {input.read1_p1} -I {input.read2_p1} -o {output.dpread1_p1} -O {output.dpread2_p1} --unpaired1 {output.duread1_p1} --unpaired2 {output.duread2_p1} &>> {log}

		## trim low quality bases
		fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2_p1} -i {output.dpread1_p1} -I {output.dpread2_p1} -o {output.pread1_p1} -O {output.pread2_p1} &>> {log}
		fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU1_p1} -i {output.duread1_p1} -o {output.uread1_p1} &>> {log}
		fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU2_p1} -i {output.duread2_p1} -o {output.uread2_p1} &>> {log}

		# population two

		## deduplicate and correct
		fastp -u {params.unqualLimit} -q {params.qualThresh} --dedup --correction -i {input.read1_p2} -I {input.read2_p2} -o {output.dpread1_p2} -O {output.dpread2_p2} --unpaired1 {output.duread1_p2} --unpaired2 {output.duread2_p2} &>> {log}

		## trim low quality bases
		fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2_p2} -i {output.dpread1_p2} -I {output.dpread2_p2} -o {output.pread1_p2} -O {output.pread2_p2} &>> {log}
		fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU1_p2} -i {output.duread1_p2} -o {output.uread1_p2} &>> {log}
		fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU2_p2} -i {output.duread2_p2} -o {output.uread2_p2} &>> {log}
		"""
