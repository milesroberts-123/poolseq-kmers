rule bwa_mem:
    input:
        amb="seqkit_results/ref_{SID}.fasta.amb",
        ann="seqkit_results/ref_{SID}.fasta.ann",
        bwt="seqkit_results/ref_{SID}.fasta.bwt",
        pac="seqkit_results/ref_{SID}.fasta.pac",
        sa="seqkit_results/ref_{SID}.fasta.sa",
        reffasta="seqkit_results/ref_{SID}.fasta",
        read1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        read2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
    output:
        bam=temp("bwa_results/{SID}_{PID}.bam"),
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_mem/{SID}_{PID}.log",
    benchmark:
        "benchmarks/bwa_mem/{SID}_{PID}.bench"
    shell:
        """
        # align reads to reference
        bwa mem -R '@RG\\tID:{wildcards.SID}\\tSM:{wildcards.SID}' -t {threads} {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output.bam}
        """
