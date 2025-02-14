rule fastp_one_pop:
	input:
		read1 = "reads_{ID}_R1.fastq",
		read2 = "reads_{ID}_R2.fastq"
	output:
		dpread1 = temp("dedup_paired_R1_{ID}.fastq"),
		dpread2 = temp("dedup_paired_R2_{ID}.fastq"),
		duread1 = temp("dedup_unpaired_R1_{ID}.fastq"),
		duread2 = temp("dedup_unpaired_R2_{ID}.fastq"),
		pread1 = temp("trimmed_paired_R1_{ID}.fastq"),
		pread2 = temp("trimmed_paired_R2_{ID}.fastq"),
		uread1 = temp("trimmed_unpaired_R1_{ID}.fastq"),
		uread2 = temp("trimmed_unpaired_R2_{ID}.fastq"),
		json = "{ID}.json"
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
		# remove duplicates, do read correction, drop low quality reads
		fastp -u 40 -q 25 --dedup --correction --json {output.json} -i {input.read1} -I {input.read2} -o {output.dpread1} -O {output.dpread2} --unpaired1 {output.duread1} --unpaired2 {output.duread2} &> {log}

		# trim low quality bases
		fastp -Q -l 31 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 30 -i {output.dpread1} -I {output.dpread2} -o {output.pread1} -O {output.pread2}
		fastp -Q -l 31 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 30 -i {output.duread1} -o {output.uread1}
		fastp -Q -l 31 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 30 -i {output.duread2} -o {output.uread2}
		"""
