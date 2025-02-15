rule hetmers_one_pop:
	input:
		"kmer_counts_{ID}.txt"
	output:
		"hetmers_{ID}_seqs.csv",
		"hetmers_{ID}_counts.csv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/hetmers.yaml"
	log: 
		"logs/hetmers/{ID}.log"
	benchmark:
		"benchmarks/hetmers/{ID}.bench"
	shell:
		"python hetmers.py -i {input} -a 2 -m 5 -o hetmers_{wildcards.ID} &> {log}"