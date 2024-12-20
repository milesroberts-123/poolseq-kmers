rule test:
	input:
		"../config/parameters.tsv"
	output:
		"slim_{ID}.vcf"
	shell:
		"touch {output}"
