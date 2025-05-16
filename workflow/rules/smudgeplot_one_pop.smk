rule smudgeplot_one_pop:
	input:
		"kmc_results/kmer_counts_{ID}.txt"
	output:
		"smudgeplot_results/{ID}_coverages.tsv",
		"smudgeplot_results/{ID}_sequences.tsv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/smudgeplot.yaml"
	log: 
		"logs/smudgeplot/{ID}.log"
	benchmark:
		"benchmarks/smudgeplot/{ID}.bench"
	shell:
		"smudgeplot.py hetkmers -o smudgeplot_results/{wildcards.ID} --middle {input} &> {log}"
