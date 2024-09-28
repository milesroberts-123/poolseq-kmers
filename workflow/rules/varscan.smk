rule varscan:
	input:
		reffasta = "ref_{ID}.fasta",
		trimbam = "trimmed.bam"
	output:
		"calls_{ID}.tsv"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp > {output}
		"""
