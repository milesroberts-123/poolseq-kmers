rule unitig_caller:
	input:
		"kmer_counts_{ID}.txt"
	output:
		unitigs_fasta = temp("unitigs_{ID}.fasta"),
		readfile = temp("reads_for_unitig-caller_{ID}.txt"),
		unitigs_rtab = temp("unitigs_{ID}.rtab"),
		tmp_fasta = temp("kmer_seqs_{ID}.fa")
	conda:
		"../envs/unitig-caller.yaml"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	params:
		rtab_prefix = "unitigs_{ID}"
	benchmark:
		"./benchmarks/unitig_caller/{ID}.bench"
	shell:
		r"""
		# turn k-mer counts into fasta
		cut -f 1 {input} | sed 's/^/>foobar\n/g' > {output.tmp_fasta}

		# create list of reads for unitig caller
		echo $PWD/{output.tmp_fasta} >> {output.readfile}

		# call unitigs
		unitig-caller --call --refs {output.readfile} --rtab --out {params.rtab_prefix}

		# convert tab output to fasta-like format
		# get just first column, remove header, add "foo" id and newline to beginning of each sequence 
		cut -f 1 {output.unitigs_rtab} | tail -n +2 | sed 's:^:>foo\n:g' > {output.unitigs_fasta}
		"""