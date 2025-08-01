rule bwa_index:
    input:
        reffasta="seqkit_results/ref_{ID}.fasta",
    output:
        amb=temp("seqkit_results/ref_{ID}.fasta.amb"),
        ann=temp("seqkit_results/ref_{ID}.fasta.ann"),
        bwt=temp("seqkit_results/ref_{ID}.fasta.bwt"),
        pac=temp("seqkit_results/ref_{ID}.fasta.pac"),
        sa=temp("seqkit_results/ref_{ID}.fasta.sa"),
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index/{ID}.log",
    benchmark:
        "benchmarks/bwa_index/{ID}.bench"
    shell:
        """
        # index reference
        bwa index {input.reffasta} &>> {log}
        """
