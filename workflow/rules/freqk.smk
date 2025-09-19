rule freqk_index:
    input:
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        fasta="ancestral_genome_results/{SID}.fasta"
    output:
        "freqk_indices/{SID}_{PID}.txt"
    benchmark:
        "benchmarks/freqk_index/{SID}_{PID}.bench"
    params:
        k=config["k"]
    shell:
        """
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output}
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
        counts="freqk_results/{SID}_{PID}_counts.txt",
        freqs="freqk_results/{SID}_{PID}_freqs.txt"
    benchmark:
        "benchmarks/freqk_count/{SID}_{PID}.bench"
    shell:
        """
        ./scripts/freqk count --index {input.index} --reads {input.pread1},{input.pread2},{input.uread1},{input.uread2} --freq-output {output.freqs} --count-output {output.counts}
        """
