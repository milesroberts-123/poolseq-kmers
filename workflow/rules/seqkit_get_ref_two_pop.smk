rule seqkit_get_ref_two_pop:
    input:
        slimfasta = "slim_results/{ID}.fasta"
    output:
        tempsamplefasta = temp("seqkit_results/samples_{ID}_p1p2.fasta"),
        p1 = "seqkit_results/samples_{ID}_p1.fasta",
        p2 = "seqkit_results/samples_{ID}_p2.fasta",
        reffasta = "seqkit_results/ref_{ID}_p1.fasta",
    params:
        n=lookup(query="ID == '{ID}'", within=parameters, cols="n")
    conda:
        "../envs/seqkit.yaml"
    log: 
        "logs/seqkit/{ID}.log"
    shell:
        """
        # create file with reference genome removed
        seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.tempsamplefasta}

        # get a haploid reference genome
        seqkit grep -n -p 1 {input.slimfasta} > {output.reffasta}

        # get first n individuals (population 1)
        seqkit head -n {params.n} {output.tempsamplefasta} > {output.p1}

        # get last n individuals (population 2)
        seqkit range -r -{params.n}:-1 {output.tempsamplefasta} > {output.p2}        
        """
