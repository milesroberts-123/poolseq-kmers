rule seqkit_get_ref_one_pop:
    input:
        slimfasta = "slim_results/{ID}.fasta"
    output:
        samplefasta = temp("seqkit_results/samples_{ID}.fasta"),
        reffasta = temp("seqkit_results/ref_{ID}.fasta"),
    conda:
        "../envs/seqkit.yaml"
    log: 
        "logs/seqkit_get_ref/{ID}.log"
    shell:
        """
        # create file with reference genome removed
        seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.samplefasta}

        # get a haploid reference genome
        seqkit grep -n -p 1 {input.slimfasta} > {output.reffasta}        
        """
