rule bwa_two_pop:
    input:
        reffasta = "seqkit_results/ref_{ID}_p1.fasta",
        read1_p1 = "fastp_results/trimmed_paired_R1_{ID}_p1.fastq",
        read2_p1 = "fastp_results/trimmed_paired_R2_{ID}_p1.fastq",
        read1_p2 = "fastp_results/trimmed_paired_R1_{ID}_p2.fastq",
        read2_p2 = "fastp_results/trimmed_paired_R2_{ID}_p2.fastq"
    output:
        bam_p1 = temp("bwa_results/{ID}_p1.bam"),
        bam_p2 = temp("bwa_results/{ID}_p2.bam"),
        amb = temp("seqkit_results/ref_{ID}_p1.fasta.amb"),
        ann = temp("seqkit_results/ref_{ID}_p1.fasta.ann"),
        bwt = temp("seqkit_results/ref_{ID}_p1.fasta.bwt"),
        pac = temp("seqkit_results/ref_{ID}_p1.fasta.pac"),
        sa = temp("seqkit_results/ref_{ID}_p1.fasta.sa")
    threads: 2
    resources:
        mem_mb_per_cpu=8000,
        time=239
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/bwa_full/{ID}.log"
    shell:
        """
        # index reference
        bwa index {input.reffasta} &>> {log}

        # align reads to reference
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1_p1} {input.read2_p1} | samtools sort -O bam > {output.bam_p1}

        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1_p2} {input.read2_p2} | samtools sort -O bam > {output.bam_p2}
        """
