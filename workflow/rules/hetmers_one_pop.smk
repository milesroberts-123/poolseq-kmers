rule hetmers_one_pop:
	input:
		"kmc_results/kmer_counts_{ID}.txt"
	output:
		"hetmers_results/{ID}_counts.csv",
		"hetmers_results/{ID}_empirical_freqs.csv",
		"hetmers_results/{ID}_bayes_states.csv",
		"hetmers_results/{ID}_hashes.csv",
		"hetmers_results/{ID}_seqs.csv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	log: 
		"logs/hetmers/{ID}.log"
	benchmark:
		"benchmarks/hetmers/{ID}.bench"
	params:
		mincount = config["mincount"],
		pool = get_pool,
		cov = get_cov
	shell:
		"""
		# create directory
		if [ ! -d "hetmers_results" ]; then
			mkdir hetmers_results
		fi

		./scripts/hetmers --inputs {input} --outputs {wildcards.ID} --coverages {params.cov} --pools {params.pool} --alphas 1 --betas 1 --minimums {params.mincount} &> {log}

		# move output to directory
		mv {wildcards.ID}_counts.csv hetmers_results/
		mv {wildcards.ID}_empirical_freqs.csv hetmers_results/
		mv {wildcards.ID}_bayes_states.csv hetmers_results/
		mv {wildcards.ID}_hashes.csv hetmers_results/
		mv {wildcards.ID}_seqs.csv hetmers_results/
                """
