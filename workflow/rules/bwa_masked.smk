rule bwa_masked:
    input:
        reffasta = "ref_masked_{ID}.fasta",
        read1 = "trimmed_paired_R1_{ID}.fastq",
        read2 = "trimmed_paired_R2_{ID}.fastq"
    output:
        bam = temp("trimmed_masked_{ID}.bam"),
        amb = temp("ref_masked_{ID}.fasta.amb"),
        ann = temp("ref_masked_{ID}.fasta.ann"),
        bwt = temp("ref_masked_{ID}.fasta.bwt"),
        pac = temp("ref_masked_{ID}.fasta.pac"),
        sa = temp("ref_masked_{ID}.fasta.sa")
    threads: 2
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/bwa_masked/{ID}.log"
    benchmark:
        "benchmarks/bwa_masked/{ID}.bench"
    shell:
        """
        # index reference
        bwa index {input.reffasta} &>> {log}

        # align reads to reference
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output.bam}
        """
