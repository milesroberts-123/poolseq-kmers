rule ancestral_samtools_faidx:
    input:
        "ancestral_genome_results/{SID}.fasta"
    output:
        temp("ancestral_genome_results/{SID}.fasta.fai")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # index reference genome
        samtools faidx {input}
        """

rule freqk_index:
    input:
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        fasta="ancestral_genome_results/{SID}.fasta",
        fai="ancestral_genome_results/{SID}.fasta.fai"
    output:
        index="freqk_indices/{SID}_{PID}.txt"
    benchmark:
        "benchmarks/freqk_index/{SID}_{PID}.bench"
    params:
        k=config["k"]
    shell:
        """
        # index panel of variants
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """

rule freqk_dedup:
    input:
        "freqk_indices/{SID}_{PID}.txt"
    output:
        "freqk_dedup/{SID}_{PID}.txt"
    benchmark:
        "benchmarks/freqk_dedup/{SID}_{PID}.bench"
    shell:
        """
        ./scripts/freqk dedup --index {input} --output {output}
        """

rule freqk_count:
    input:
        index="freqk_dedup/{SID}_{PID}.txt",
        pread1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
        uread1="fastp_results/trimmed_unpaired_R1_{SID}_{PID}.fastq",
        uread2="fastp_results/trimmed_unpaired_R2_{SID}_{PID}.fastq",
    output:
        all=temp("all_{SID}_{PID}.fastq"),
        counts="freqk_results/{SID}_{PID}_counts.txt",
        freqs="freqk_results/{SID}_{PID}_freqs.txt"
    benchmark:
        "benchmarks/freqk_count/{SID}_{PID}.bench"
    shell:
        """
        cat {input.pread1} {input.pread2} {input.uread1} {input.uread2} > {output.all}
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {output.all} --freq-output {output.freqs} --count-output {output.counts}
        """
