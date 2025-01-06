rule bedtools:
	input:
		fasta = "ref_{ID}.fasta",
		bed = "../config/mask.bed"
	output:
		"ref_masked_{ID}.fasta"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/bedtools.yaml"
	log: 
		"logs/bedtools/{ID}.log"
	shell:
		"""
		bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}
		"""