rule bwa_one_pop:
    input:
        reffasta="seqkit_results/ref_{ID}.fasta",
        read1="fastp_results/trimmed_paired_R1_{ID}.fastq",
        read2="fastp_results/trimmed_paired_R2_{ID}.fastq",
    output:
        bam=temp("bwa_results/{ID}.bam"),
        amb=temp("seqkit_results/ref_{ID}.fasta.amb"),
        ann=temp("seqkit_results/ref_{ID}.fasta.ann"),
        bwt=temp("seqkit_results/ref_{ID}.fasta.bwt"),
        pac=temp("seqkit_results/ref_{ID}.fasta.pac"),
        sa=temp("seqkit_results/ref_{ID}.fasta.sa"),
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_full/{ID}.log",
    benchmark:
        "benchmarks/bwa/{ID}.bench"
    shell:
        """
        # index reference
        bwa index {input.reffasta} &>> {log}

        # align reads to reference
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output.bam}
        """
