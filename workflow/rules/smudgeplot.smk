rule smudgeplot:
	input:
		"kmer_counts_{ID}.txt"
	output:
		"kmerpairs_{ID}_coverages.tsv",
		"kmerpairs_{ID}_sequences.tsv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/smudgeplot.yaml"
	log: 
		"logs/smudgeplot/{ID}.log"
	shell:
		"smudgeplot.py hetkmers -o kmerpairs_{wildcards.ID} --middle {input} &> {log}"
