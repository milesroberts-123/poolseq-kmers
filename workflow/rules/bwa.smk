rule bwa:
	input:
		reffasta = "ref_{ID}.fasta",
		read1 = "trimmed_paired_R1_{ID}.fastq",
		read2 = "trimmed_paired_R2_{ID}.fastq"
	output:
		"trimmed_{ID}.bam"
	threads: 2
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bwa.yaml"
	log: 
		"logs/bwa/{ID}.log"
	shell:
		"""
		# index reference
		bwa index {input.reffasta} &>> {log}

		# align reads to reference
		bwa mem -R '@RG\tID:{wildcards.ID}\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output}
		"""
