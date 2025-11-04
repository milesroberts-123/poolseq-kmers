rule ref_samtools_faidx:
    input:
        ref="seqkit_results/ref_{ID}.fasta",
    output:
        fai=temp("seqkit_results/ref_{ID}.fasta.fai"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # index reference
        samtools faidx {input.ref}
        """
