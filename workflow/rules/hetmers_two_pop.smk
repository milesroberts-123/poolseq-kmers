rule hetmers_two_pop:
	input:
		p1="kmer_counts_{ID}_p1.txt",
		p2="kmer_counts_{ID}_p2.txt"
	output:
		"hetmers_{ID}_p1_counts.csv",
		"hetmers_{ID}_p1_empirical_freqs.csv",
		"hetmers_{ID}_p1_hashes.csv",
		"hetmers_{ID}_p1_seqs.csv",
		"hetmers_{ID}_p2_counts.csv",
		"hetmers_{ID}_p2_empirical_freqs.csv",
		"hetmers_{ID}_p2_hashes.csv",
		"hetmers_{ID}_p2_seqs.csv"
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
		./scripts/hetmers --inputs {input.p1} --outputs hetmers_{wildcards.ID}_p1 --analyses hetmers --coverages 1 --pools 1 --alphas 1 --betas 1 --minimums {params.mincount} &> {log}

		./scripts/hetmers --inputs {input.p2} --outputs hetmers_{wildcards.ID}_p2 --analyses hetmers --coverages 1 --pools 1 --alphas 1 --betas 1 --minimums {params.mincount} &> {log}
                """
