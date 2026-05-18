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

rule seqkit_rename:
    input:
        "unitig_caller_results/unitigs_{SID}_{PID}.fasta",
    output:
        "unitig_caller_results/unitigs_renamed_{SID}_{PID}.fasta",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit rename {input} > {output}
        """

rule makeblastdb:
    input:
        ref="ancestral_genome_results/{SID}.fasta",
    output:
        ndb=temp("{SID}.ndb"),
        nhr=temp("{SID}.nhr"),
        nin=temp("{SID}.nin"),
        njs=temp("{SID}.njs"),
        not_dbfile=temp("{SID}.not"),
        nsq=temp("{SID}.nsq"),
        ntf=temp("{SID}.ntf"),
        nto=temp("{SID}.nto"),
    conda:
        "../envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input.ref} -title $(basename {input.ref} .fasta) -dbtype nucl -out $(basename {input.ref} .fasta)
        """

rule blast:
    input:
        ref="ancestral_genome_results/{SID}.fasta",
        ndb="{SID}.ndb",
        nhr="{SID}.nhr",
        nin="{SID}.nin",
        njs="{SID}.njs",
        not_dbfile="{SID}.not",
        nsq="{SID}.nsq",
        ntf="{SID}.ntf",
        nto="{SID}.nto",
        unitigs="unitig_caller_results/unitigs_renamed_{SID}_{PID}.fasta",
    output:
        alignments="blast_results/{SID}_{PID}.txt",
    conda:
        "../envs/blast.yaml"
    params:
        blastEvalue=config["blastEvalue"],
    shell:
        """
        blastn -query {input.unitigs} -db $(basename {input.ref} .fasta) -out {output.alignments} -max_target_seqs 1 -evalue {params.blastEvalue} -outfmt 6
        """
