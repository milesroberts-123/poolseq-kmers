rule seqkit_rename:
	input:
		"unitigs_{ID}.fasta"
	output:
		"unitigs_renamed_{ID}.fasta"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/seqkit.yaml"	
	shell:
		"""
		seqkit rename {input} > {output}
		"""