rule varscan:
	input:
		reffasta = "ref_{ID}.fasta",
		trimbam = "trimmed_{ID}.bam"
	output:
		"calls_{ID}.tsv"
	threads: 8
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/varscan.yml"
	log: 
		"logs/varscan_{ID}.log"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp > {output} &> {log}
		"""
