rule seqkit_get_ref:
    input:
        slimfasta="slim_results/{ID}.fasta",
    output:
        reffasta="seqkit_results/ref_{ID}.fasta",
        samplefasta="seqkit_results/samples_{ID}.fasta",
    conda:
        "../envs/seqkit.yaml"
    log:
        "logs/seqkit/{ID}.log",
    shell:
        """
        # get a haploid reference genome
        seqkit grep -n -p 1 {input.slimfasta} > {output.reffasta}

        # create file with reference genome removed
        seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.samplefasta}
        """
