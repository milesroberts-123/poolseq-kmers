rule varscan_full:
	input:
		reffasta = "ref_full_{ID}.fasta",
		trimbam = "trimmed_full_{ID}.bam"
	output:
		"calls_full_{ID}.tsv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/varscan.yaml"
	log: 
		"logs/varscan_full/{ID}.log"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp 1> {output} 2> {log}
		"""
