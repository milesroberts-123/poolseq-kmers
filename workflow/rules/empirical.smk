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
        r1=expand("real_fastp_results/trimmed_paired_R1_{ID}.fastq", ID=config["real_accessions"]),
        r2=expand("real_fastp_results/trimmed_paired_R2_{ID}.fastq", ID=config["real_accessions"]),
        u1=expand("real_fastp_results/trimmed_unpaired_R1_{ID}.fastq", ID=config["real_accessions"]),
        u2=expand("real_fastp_results/trimmed_unpaired_R2_{ID}.fastq", ID=config["real_accessions"]),
    output:
        all="pseudopool.fastq",
        r1="r1_pool.fastq",
        r2="r2_pool.fastq",
        u="u_pool.fastq",
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        cat {input.u1} {input.u2} > {output.u}
        cat {output.r1} {output.r2} {output.u} > {output.all}
        """


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
        ./scripts/freqk ref-dedup --index {input.index} --fasta {input.fasta} --vcf {input.vcf} --output {output}
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

rule real_vg_autoindex:
    input:
        vcf=config["real_vcf"],
        fasta=config["real_fasta"],
    output:
        dist = temp("real.dist"),
        gbz = temp("real.giraffe.gbz"),
        min = temp("real.shortread.withzip.min"),
        zip = temp("real.shortread.zipcodes")
    benchmark:
        "benchmarks/real_data/vg_autoindex.bench"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe -r {input.fasta} -v {input.vcf} -p real"

rule real_vg_giraffe_paired:
    input:
        dist = "real.dist",
        gbz = "real.giraffe.gbz",
        min = "real.shortread.withzip.min",
        zip = "real.shortread.zipcodes",
        pread1 = "r1_pool.fastq",
        pread2 = "r2_pool.fastq",
    output:
        temp("real_vg_giraffe_results/paired.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/real_data/vg_giraffe_paired.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -p -t {threads} --rescue-algorithm none -f {input.pread1} -f {input.pread2} > {output}
        """

rule real_vg_giraffe_unpaired:
    input:
        dist = "real.dist",
        gbz = "real.giraffe.gbz",
        min = "real.shortread.withzip.min",
        zip = "real.shortread.zipcodes",
        uread = "u_pool.fastq",
    output:
        temp("real_vg_giraffe_results/unpaired.gam")
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/real_data/vg_giraffe_unpaired.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} --rescue-algorithm none -p -t {threads} -f {input.uread} > {output}
        """

rule real_vg_surject:
    input:
        gbz = "real.giraffe.gbz",
        gam = "real_vg_giraffe_results/{pairing}.gam",
    output:
        temp("real_vg_surject_results/{pairing}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/real_data/vg_surject_{pairing}.bench"
    shell:
        """
        vg surject -x {input.gbz} --progress -t {threads} -b {input.gam} > {output}
        """

rule real_samtools_sort:
    input:
        "real_vg_surject_results/{pairing}.bam"
    output:
        temp("real_samtools_sort_results/{pairing}.bam")
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/real_data/samtools_sort_{pairing}.bench"
    shell:
        """
        samtools sort {input} -o {output}
        """

rule real_samtools_merge:
    input:
        expand("real_samtools_sort_results/{pairing}.bam", pairing = ["unpaired", "paired"])
    output:
        temp("merged.bam")
    benchmark:
        "benchmarks/real_data/samtools_merge.bench"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools merge -o {output} {input}
        """

rule split_merged_bam:
    input:
        "merged.bam"
    output:
        "splits/{chr}.bam"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "samtools view -b {input} {wildcards.chr} > {output}"

rule real_freebayes_vg:
    input:
        reffasta=config["real_fasta"],
        trimbam="splits/{chr}.bam",
        #vcf=config["real_vcf"],
        vcf="{chr}.vcf.gz"
    output:
        "calls_{chr}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/real_data/freebayes_vg_{chr}.bench"
    params:
        n=config["poolsize"]
    shell:
        """
        freebayes -f {input.reffasta} -p {params.n} --min-alternate-count 2 --min-alternate-fraction 0.001 --use-best-n-alleles 4 -g 1000 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} > {output}
        """

