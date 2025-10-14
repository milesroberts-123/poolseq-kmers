rule downsample_reads:
    input:
        r1="../config/real_reads/{sample}_R1.fastq.gz",
        r2="../config/real_reads/{sample}_R2.fastq.gz"
    output:
        r1="downsampled_reads/{sample}_R1.fastq.gz",
        r2="downsampled_reads/{sample}_R2.fastq.gz"
    conda:
        "../envs/seqtk.yaml"
    params:
        downfactor
    shell:
        """
        """

rule pseudopool:

rule real_fastp:
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        """

rule real_freqk_index:
    input:
        vcf=config["real_vcf"],
        fasta=config["real_fasta"],
    output:
        index="freqk_indices/{species}.txt"
    benchmark:
        "benchmarks/freqk_index/{species}.bench"
    params:
        k=config["k"]
    shell:
        """
        # index panel of variants
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """

rule real_freqk_dedup:
    input:
        "freqk_indices/{species}.txt"
    output:
        "freqk_dedup/{species}.txt"
    benchmark:
        "benchmarks/freqk_dedup/{species}.bench"
    shell:
        """
        ./scripts/freqk dedup --index {input} --output {output}
        """

rule real_freqk_count:
    input:
        index="freqk_dedup/{species}.txt",
        pread1="fastp_results/trimmed_paired_R1_{species}.fastq.gz",
        pread2="fastp_results/trimmed_paired_R2_{species}.fastq.gz",
        uread1="fastp_results/trimmed_unpaired_R1_{species}.fastq.gz",
        uread2="fastp_results/trimmed_unpaired_R2_{species}.fastq.gz",
    output:
        all=temp("all_{species}.fastq"),
        counts="freqk_results/{species}_counts.txt",
        freqs="freqk_results/{species}_freqs.txt"
    benchmark:
        "benchmarks/freqk_count/{species}.bench"
    shell:
        """
        cat {input.pread1} {input.pread2} {input.uread1} {input.uread2} > {output.all}
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {output.all} --freq-output {output.freqs} --count-output {output.counts}
        """


