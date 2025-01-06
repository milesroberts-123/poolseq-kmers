rule unitig_caller:
	input:
		pread1 = "trimmed_paired_R1_{ID}.fastq",
		pread2 = "trimmed_paired_R2_{ID}.fastq",
		uread1 = "trimmed_unpaired_R1_{ID}.fastq",
		uread2 = "trimmed_unpaired_R2_{ID}.fastq"
	output:
		unitigs_fasta = temp("unitigs_{ID}.fasta"),
		readfile = temp("reads_for_unitig-caller_{ID}.txt"),
		unitigs_rtab = temp("unitigs_{ID}.rtab")
	conda:
		"../envs/unitig-caller.yaml"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	params:
		rtab_prefix = "unitigs_{ID}"
	shell:
		r"""
		# create list of reads for unitig caller
		echo $PWD/{input.pread1} >> {output.readfile}
		echo $PWD/{input.pread2} >> {output.readfile}
		echo $PWD/{input.uread1} >> {output.readfile}
		echo $PWD/{input.uread2} >> {output.readfile}

		# call unitigs
		unitig-caller --call --reads {output.readfile} --rtab --out {params.rtab_prefix}

		# convert tab output to fasta-like format
		# get just first column, remove header, add "foo" id and newline to beginning of each sequence 
		cut -f 1 {output.unitigs_rtab} | tail -n +2 | sed 's:^:>foo\n:g' > {output.unitigs_fasta}
		"""