rule sra:
    output:
        r1="{ID}_1.fastq",
        r2="{ID}_2.fastq"
    retries: 2
    conda:
        "../envs/sra.yaml"
    params:
        k=config["k"],
        N=config["num_reads"],
    shell:
        """
        fastq-dump --split-3 --skip-technical -X {params.N} --clip -M {params.k} {wildcards.ID}
        #fasterq-dump --threads {threads} --split-files --skip-technical {wildcards.ID}
        """

rule real_fastp:
    input:
        r1="{ID}_1.fastq",
        r2="{ID}_2.fastq"    
    output:
        pread1=temp("real_fastp_results/trimmed_paired_R1_{ID}.fastq"),
        pread2=temp("real_fastp_results/trimmed_paired_R2_{ID}.fastq"),
        uread1=temp("real_fastp_results/trimmed_unpaired_R1_{ID}.fastq"),
        uread2=temp("real_fastp_results/trimmed_unpaired_R2_{ID}.fastq"),
        jsonR1R2="real_fastp_results/{ID}_R1R2.json",
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
    shell:
        """
        fastp --thread {threads} -u {params.unqualLimit} -q {params.qualThresh} --correction -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {input.r1} -I {input.r2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2}
        """

#rule downsample:
#    input:
#        "real_fastp_results/trimmed_paired_R1_{ID}.fastq",
#        "real_fastp_results/trimmed_paired_R2_{ID}.fastq",
#        "real_fastp_results/trimmed_unpaired_R1_{ID}.fastq",
#        "real_fastp_results/trimmed_unpaired_R2_{ID}.fastq",
#    output:
#        temp("downsample/{ID}.fastq")
#    conda:
#        "../envs/seqkit.yaml"
#    params:
#        N=config["num_reads"]
#    shell:
#        """
#        set +o pipefail; cat {input} | seqkit sample -s 21 -p 0.1 | seqkit head -n {params.N} > {output}
#        """

rule real_cat:
    input:
        #expand("downsample/{ID}.fastq", ID=config["real_accessions"]),
        expand("real_fastp_results/trimmed_paired_R1_{ID}.fastq", ID=config["real_accessions"]),
        expand("real_fastp_results/trimmed_paired_R2_{ID}.fastq", ID=config["real_accessions"]),
        expand("real_fastp_results/trimmed_unpaired_R1_{ID}.fastq", ID=config["real_accessions"]),
        expand("real_fastp_results/trimmed_unpaired_R2_{ID}.fastq", ID=config["real_accessions"]),
    output:
        "pseudopool.fastq"
    shell:
        "cat {input} > {output}"


rule real_freqk_index:
    input:
        vcf=config["real_vcf"],
        fasta=config["real_fasta"],
    output:
        index=temp("index.txt")
    benchmark:
        "benchmarks/real_data/freqk_index.bench"
    params:
        k=config["k"]
    shell:
        """
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """

rule real_freqk_var_dedup:
    input:
        "index.txt"
    output:
        temp("var_index.txt")
    benchmark:
        "benchmarks/real_data/freqk_var_dedup.bench"
    shell:
        """
        ./scripts/freqk var-dedup --index {input} --output {output}
        """

rule real_freqk_ref_dedup:
    input:
        index="var_index.txt",
        vcf=config["real_vcf"],
        fasta=config["real_fasta"],
    output:
        "ref_index.txt"
    benchmark:
        "benchmarks/real_data/freqk_ref_dedup.bench"
    params:
        k=config["k"]
    shell:
        """
        ./scripts/freqk ref-dedup --index {input.index} --fasta {input.fasta} --vcf {input.vcf} --kmer {params.k} --output {output}
        """

rule real_freqk_count:
    input:
        reads="pseudopool.fastq",
        index="ref_index.txt",
    output:
        counts="real_freqk_results/counts_by_kmer.txt",
        freqs="real_freqk_results/counts_by_allele.txt"
    benchmark:
        "benchmarks/real_data/freqk_count.bench"
    shell:
        """
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {input.reads} --freq-output {output.freqs} --count-output {output.counts}
        """

rule real_freqk_call:
    input:
        counts="real_freqk_results/counts_by_allele.txt",
        index="ref_index.txt"
    output:
        "calls.txt"
    benchmark:
        "benchmarks/real_data/freqk_call.bench"
    shell:
        "./scripts/freqk call --index {input.index} -c {input.counts} --output {output}"
