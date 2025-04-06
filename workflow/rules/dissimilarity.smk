rule dissimilarity:
	input:
		p1="kmer_counts_{ID}_p1.txt",
		p2="kmer_counts_{ID}_p2.txt"
	output:
		"dissimilarity_{ID}.txt"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	log: 
		"logs/dissimilarity/{ID}.log"
	conda:
		"../envs/R.yaml"
	shell:
		"""
		Rscript scripts/dissimilarity.R {input.p1} {input.p2} {output} {threads} &> {log}
                """
