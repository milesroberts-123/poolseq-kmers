rule smudgeplot:
	input:
		"kmer_counts_{ID}.txt"
	output:
	
	shell:
		"smudgeplot.py hetkmers --middle {input}"
