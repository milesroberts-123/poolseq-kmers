rule ref_samtools_faidx:
    input:
        ref="seqkit_results/ref_{ID}.fasta",
    output:
        fai=temp("seqkit_results/ref_{ID}.fasta.fai"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/samtools_faidx/{ID}.log",
    shell:
        """
        # index reference
        samtools faidx {input.ref}
        """
