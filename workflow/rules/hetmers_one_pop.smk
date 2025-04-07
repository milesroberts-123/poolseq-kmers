rule hetmers_one_pop:
	input:
		"kmer_counts_{ID}.txt"
	output:
		"hetmers_{ID}_counts.csv",
		"hetmers_{ID}_empirical_freqs.csv",
		"hetmers_{ID}_hashes.csv",
		"hetmers_{ID}_seqs.csv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	log: 
		"logs/hetmers/{ID}.log"
	benchmark:
		"benchmarks/hetmers/{ID}.bench"
	params:
		mincount = config["mincount"]
	shell:
		"""
		./scripts/hetmers --inputs {input} --outputs hetmers_{wildcards.ID} --coverages 1 --pools 1 --alphas 1 --betas 1 --minimums {params.mincount} &> {log}
                """