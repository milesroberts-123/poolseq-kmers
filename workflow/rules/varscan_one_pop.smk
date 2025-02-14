rule varscan_one_pop:
	input:
		reffasta = "ref_{ID}.fasta",
		trimbam = "trimmed_{ID}.bam"
	output:
		"calls_{ID}.tsv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/varscan.yaml"
	log: 
		"logs/varscan/{ID}.log"
	benchmark:
		"benchmarks/varscan/{ID}.bench"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp 1> {output} 2> {log}
		"""
