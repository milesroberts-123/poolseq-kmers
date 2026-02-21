rule sra:
    output:
        r1="{ID}_1.fastq",
        r2="{ID}_2.fastq",
    retries: 2
    conda:
        "../envs/sra.yaml"
    params:
        N=config["num_reads"],
    shell:
        """
        fasterq-dump --threads {threads} --split-files --skip-technical {wildcards.ID}
        """


rule downsample:
    input:
        r1="{ID}_1.fastq",
        r2="{ID}_2.fastq",
    output:
        r1=temp("downsample/{ID}_{num}_1.fastq"),
        r2=temp("downsample/{ID}_{num}_2.fastq"),
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        set +o pipefail
        seqkit sample -s 21 -p 0.1 {input.r1} | seqkit head -n {wildcards.num} > {output.r1}
        seqkit sample -s 21 -p 0.1 {input.r2} | seqkit head -n {wildcards.num} > {output.r2}
        """

rule real_fastp:
    input:
        r1="downsample/{ID}_{num}_1.fastq",
        r2="downsample/{ID}_{num}_2.fastq",
    output:
        pread1=temp("real_fastp_results/trimmed_paired_R1_{ID}_{num}.fastq"),
        pread2=temp("real_fastp_results/trimmed_paired_R2_{ID}_{num}.fastq"),
        uread1=temp("real_fastp_results/trimmed_unpaired_R1_{ID}_{num}.fastq"),
        uread2=temp("real_fastp_results/trimmed_unpaired_R2_{ID}_{num}.fastq"),
        jsonR1R2="real_fastp_results/{ID}_{num}_R1R2.json",
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=config["real_k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
    shell:
        """
        fastp --thread {threads} -u {params.unqualLimit} -q {params.qualThresh} --correction -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {input.r1} -I {input.r2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2}
        """


rule real_cat:
    input:
        r1=expand(
            "real_fastp_results/trimmed_paired_R1_{ID}_{{num}}.fastq",
            ID=config["real_accessions"],
        ),
        r2=expand(
            "real_fastp_results/trimmed_paired_R2_{ID}_{{num}}.fastq",
            ID=config["real_accessions"],
        ),
        u1=expand(
            "real_fastp_results/trimmed_unpaired_R1_{ID}_{{num}}.fastq",
            ID=config["real_accessions"],
        ),
        u2=expand(
            "real_fastp_results/trimmed_unpaired_R2_{ID}_{{num}}.fastq",
            ID=config["real_accessions"],
        ),
    output:
        all=temp("pseudopool_{num}.fastq"),
        r1=temp("r1_pool_{num}.fastq"),
        r2=temp("r2_pool_{num}.fastq"),
        u=temp("u_pool_{num}.fastq"),
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        cat {input.u1} {input.u2} > {output.u}
        cat {output.r1} {output.r2} {output.u} > {output.all}
        """


rule real_samtools_faidx:
    input:
        config["real_fasta"],
    output:
        config["real_fasta"] + ".fai",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools faidx {input}
        """


rule real_freqk_index:
    input:
        vcf=config["real_vcf"],
        fasta=config["real_fasta"],
    output:
        index=temp("index.txt"),
    benchmark:
        "benchmarks/real_data/freqk_index.bench"
    params:
        k=config["real_k"],
    shell:
        """
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """


rule real_freqk_var_dedup:
    input:
        "index.txt",
    output:
        temp("var_index.txt"),
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
        "ref_index.txt",
    benchmark:
        "benchmarks/real_data/freqk_ref_dedup.bench"
    shell:
        """
        ./scripts/freqk ref-dedup --index {input.index} --fasta {input.fasta} --vcf {input.vcf} --output {output}
        """


rule real_freqk_count:
    input:
        reads="pseudopool_{num}.fastq",
        index="ref_index.txt",
    output:
        counts="real_freqk_results/{num}/counts_by_kmer.txt",
        freqs="real_freqk_results/{num}/counts_by_allele.txt",
    benchmark:
        "benchmarks/real_data/freqk_count_{num}.bench"
    shell:
        """
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {input.reads} --freq-output {output.freqs} --count-output {output.counts}
        """


rule real_freqk_call:
    input:
        counts="real_freqk_results/{num}/counts_by_allele.txt",
        index="ref_index.txt",
    output:
        "real_freqk_results/{num}/calls.txt",
    benchmark:
        "benchmarks/real_data/freqk_call_{num}.bench"
    shell:
        "./scripts/freqk call --index {input.index} -c {input.counts} --output {output}"


rule real_vg_autoindex:
    input:
        vcf=config["real_vcf"],
        fasta=config["real_fasta"],
    output:
        dist=temp("real.dist"),
        gbz=temp("real.giraffe.gbz"),
        min=temp("real.shortread.withzip.min"),
        zip=temp("real.shortread.zipcodes"),
    benchmark:
        "benchmarks/real_data/vg_autoindex.bench"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe -r {input.fasta} -v {input.vcf} -p real"


rule real_vg_giraffe_paired:
    input:
        dist="real.dist",
        gbz="real.giraffe.gbz",
        min="real.shortread.withzip.min",
        zip="real.shortread.zipcodes",
        pread1="r1_pool_{num}.fastq",
        pread2="r2_pool_{num}.fastq",
    output:
        temp("real_vg_giraffe_results/{num}/paired.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/real_data/vg_giraffe_paired_{num}.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -p -t {threads} --rescue-algorithm none -f {input.pread1} -f {input.pread2} > {output}
        """


rule real_vg_giraffe_unpaired:
    input:
        dist="real.dist",
        gbz="real.giraffe.gbz",
        min="real.shortread.withzip.min",
        zip="real.shortread.zipcodes",
        uread="u_pool_{num}.fastq",
    output:
        temp("real_vg_giraffe_results/{num}/unpaired.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/real_data/vg_giraffe_unpaired_{num}.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} --rescue-algorithm none -p -t {threads} -f {input.uread} > {output}
        """


rule real_vg_surject:
    input:
        gbz="real.giraffe.gbz",
        gam="real_vg_giraffe_results/{num}/{pairing}.gam",
    output:
        temp("real_vg_surject_results/{num}/{pairing}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/real_data/vg_surject_{pairing}_{num}.bench"
    shell:
        """
        vg surject -x {input.gbz} --progress -t {threads} -b {input.gam} > {output}
        """


rule real_samtools_sort:
    input:
        "real_vg_surject_results/{num}/{pairing}.bam",
    output:
        temp("real_samtools_sort_results/{num}/{pairing}.bam"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/real_data/samtools_sort_{pairing}_{num}.bench"
    shell:
        """
        samtools sort {input} -o {output}
        """


rule real_samtools_merge:
    input:
        expand(
            "real_samtools_sort_results/{{num}}/{pairing}.bam",
            pairing=["unpaired", "paired"],
        ),
    output:
        temp("samtools_merge_results/{num}.bam"),
    benchmark:
        "benchmarks/real_data/samtools_merge_{num}.bench"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools merge -o {output} {input}
        """

rule real_samtools_index:
    input:
        "samtools_merge_results/{num}.bam"
    output:
        "samtools_merge_results/{num}.bam.bai"
    benchmark:
        "benchmarks/real_data/samtools_index_{num}.bench"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools index {input}
        """

rule split_merged_bam:
    input:
        "samtools_merge_results/{num}.bam",
        "samtools_merge_results/{num}.bam.bai"
    output:
        "splits/{num}/{chr}.bam",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/real_data/split_merged_bam/{num}/{chr}.bench"
    shell:
        "samtools view -b {input} {wildcards.chr} > {output}"

rule split_vcf:
    input:
        config["real_vcf"],
    output:
        "splits/{chr}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/real_data/split_vcf/{chr}.bench"
    shell:
        "bcftools view -r {wildcards.chr} -Oz -o {output} {input}"

rule index_vcf:
    input:
        "splits/{chr}.vcf.gz"
    output:
        "splits/{chr}.vcf.gz.tbi"
    benchmark:
        "benchmarks/real_data/index_vcf/{chr}.bench"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix {input}"

rule real_freebayes_vg:
    input:
        reffasta=config["real_fasta"],
        fai=config["real_fasta"] + ".fai",
        trimbam="splits/{num}/{chr}.bam",
        #vcf=config["real_vcf"],
        vcf="splits/{chr}.vcf.gz",
        tbi="splits/{chr}.vcf.gz.tbi"
    output:
        "real_freebayes_results/{num}_{chr}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/real_data/freebayes_vg_{num}_{chr}.bench"
    params:
        n=config["poolsize"],
    shell:
        """
        freebayes-parallel <(fasta_generate_regions.py {input.fai} 100000) {threads} -f {input.reffasta} -p {params.n} --min-alternate-count 2 --min-alternate-fraction 0.001 --use-best-n-alleles 4 -g 1000 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} > {output}
        #freebayes -f {input.reffasta} -p {params.n} --min-alternate-count 2 --min-alternate-fraction 0.001 --use-best-n-alleles 4 -g 1000 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} > {output}
        """

rule real_varscan_vg:
    input:
        reffasta=config["real_fasta"],
        trimbam="splits/{num}/{chr}.bam",
    output:
        "real_varscan_results/{num}_{chr}.tsv",
    conda:
        "../envs/varscan.yaml"
    benchmark:
        "benchmarks/real_data/varscan_vg_{num}_{chr}.bench"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp 1> {output}
        """

rule real_poolsnp_vg:
    input:
        reffasta=config["real_fasta"],
        trimbam="splits/{num}/{chr}.bam",
    output:
        vcf="{chr}_{num}_real_poolsnp_results.vcf.gz",
        cov=temp("{chr}_{num}_real_poolsnp_results-cov-0.9999.txt"),
        bs=temp("{chr}_{num}_real_poolsnp_results_BS.txt.gz"),
        mpileup=temp("{chr}_{num}.mpileup"),
    params:
        wd=get_wd,
        mincount=config["mincount"],
    conda:
        "../envs/poolsnp.yaml"
    benchmark:
        "benchmarks/real_data/poolsnp_vg_{num}_{chr}.bench"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} > {output.mpileup}

        PoolSNP.sh   \
        mpileup={params.wd}{output.mpileup} \
        reference={params.wd}{input.reffasta} \
        names=foobar \
        max-cov=0.9999 \
        min-cov={params.mincount} \
        min-count={params.mincount} \
        min-freq=0.01 \
        miss-frac=0 \
        badsites=1 \
        allsites=0 \
        output={params.wd}{wildcards.chr}_{wildcards.num}_real_poolsnp_results
        """
