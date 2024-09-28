rule bwa:
	input:
		reffasta = "ref_{ID}.fasta",
		read1 = "trimmed_paired_R1_{ID}.fastq",
		read2 = "trimmed_paired_R2_{ID}.fastq"
	output:
		"trimmed_{ID}.bam"
	threads: 8 
	shell:
		"""
		# index reference
		bwa index {input.reffasta}

		# align reads to reference
		bwa mem -R '@RG\tID:{wildcards.ID}\tSM:{wildcards.ID}' -t threads {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output}
		"""
