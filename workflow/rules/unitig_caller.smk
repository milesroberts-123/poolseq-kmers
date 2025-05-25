rule unitig_caller:
    input:
        "kmc_results/kmer_counts_{SID}_{PID}.txt",
    output:
        unitigs_fasta=temp("unitig_caller_results/unitigs_{SID}_{PID}.fasta"),
        readfile=temp("unitig_caller_results/reads_for_unitig-caller_{SID}_{PID}.txt"),
        unitigs_rtab=temp("unitig_caller_results/unitigs_{SID}_{PID}.rtab"),
        tmp_fasta=temp("unitig_caller_results/kmer_seqs_{SID}_{PID}.fa"),
    conda:
        "../envs/unitig-caller.yaml"
    params:
        rtab_prefix="unitig_caller_results/unitigs_{SID}_{PID}",
    benchmark:
        "benchmarks/unitig_caller/{SID}_{PID}.bench"
    log:
        "logs/unitig_caller/{SID}_{PID}.log",
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
