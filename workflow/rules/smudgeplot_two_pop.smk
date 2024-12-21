rule smudgeplot_two_pop:
	input:
		p1="kmer_counts_{ID}_p1.txt",
		p2="kmer_counts_{ID}_p2.txt"
	output:
		"kmerpairs_{ID}_p1_coverages.tsv",
		"kmerpairs_{ID}_p1_sequences.tsv",
		"kmerpairs_{ID}_p2_coverages.tsv",
		"kmerpairs_{ID}_p2_sequences.tsv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/smudgeplot.yaml"
	log: 
		"logs/smudgeplot/{ID}.log"
	shell:
		"""
		smudgeplot.py hetkmers -o kmerpairs_{wildcards.ID}_p1 --middle {input.p1} &> {log}
		smudgeplot.py hetkmers -o kmerpairs_{wildcards.ID}_p2 --middle {input.p2} &> {log}
		"""
