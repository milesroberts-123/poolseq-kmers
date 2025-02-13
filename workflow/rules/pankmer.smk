rule pankmer:
	input:
		reffasta = "ref_{ID}.fasta",
		pread1 = temp("trimmed_paired_R1_{ID}.fastq"),
		pread2 = temp("trimmed_paired_R2_{ID}.fastq"),
		uread1 = temp("trimmed_unpaired_R1_{ID}.fastq"),
		uread2 = temp("trimmed_unpaired_R2_{ID}.fastq")
	output:
		reads_index="reads_index_{ID}.tar",
		table="anchor_values_{ID}.tsv",
		plot="anchor_plot_{ID}.svg"
	conda:
		"../envs/pankmer.yaml"		
	shell:
		"""
		# create pankmer index
		pankmer index -g {input.pread1} {input.pread2} {input.uread1} {input.uread2} -o {output.reads_index}

		# create anchor plot
		pankmer anchor-genome -i {output.reads_index} \
		-a {input.reffasta} \
		-c 1 \
		--table {output.anchor_values} \
		-o {output.anchor_plot}
		"""