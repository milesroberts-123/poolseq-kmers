rule seqkit_get_ref:
    input:
        slimfasta="slim_results/{SID}.fasta",
    output:
        reffasta="seqkit_results/ref_{SID}.fasta",
        samplefasta="seqkit_results/samples_across_pop_{SID}.fasta",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # get a haploid reference genome
        seqkit grep -n -p 1 {input.slimfasta} > {output.reffasta}

        # create file with reference genome removed
        seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.samplefasta}
        """
