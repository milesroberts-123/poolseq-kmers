rule sra:
    output:
        r1=temp("{ID}_1.fastq"),
        r2=temp("{ID}_2.fastq")
    conda:
        "../envs/sra.yaml"
    shell:
        """
        fasterq-dump --threads {threads} --split-files --skip-technical {wildcards.ID}
        """

rule pseudopool:
    input:
        r1=expand("{ID}_1.fastq", ID=),
        r2=expand("{ID}_2.fastq", ID=)
    output:
        r1="pseudopools/{species}_R1.fastq.gz",
        r2="pseudopools/{species}_R2.fastq.gz"
    conda:
        "../envs/seqtk.yaml"
    params:
        N=lookup()
    shell:
        """
        zcat {input.r1} | seqkit sample -s 21 -p 0.1 | seqkit head -n {params.N} > {output.r1}
        zcat {input.r2} | seqkit sample -s 21 -p 0.1 | seqkit head -n {params.N} > {output.r2}
        """

rule real_fastp:
    input:
        r1="pseudopools/{pool}_R1.fastq.gz",
        r2="pseudopools/{pool}_R2.fastq.gz"    
    output:
        pread1=temp("real_fastp_results/trimmed_paired_R1_{pool}.fastq"),
        pread2=temp("real_fastp_results/trimmed_paired_R2_{pool}.fastq"),
        uread1=temp("real_fastp_results/trimmed_unpaired_R1_{pool}.fastq"),
        uread2=temp("real_fastp_results/trimmed_unpaired_R2_{pool}.fastq"),
        jsonR1R2="real_fastp_results/{pool}_R1R2.json",
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
    shell:
        """
        fastp --thread {threads} -u {params.unqualLimit} -q {params.qualThresh} --correction -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {input.read1} -I {input.read2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2}
        """

rule real_samtools:
    input:
        fasta=config["real_fasta"],
    output:
        config["real_fasta"] + ".fai"
    conda:
        ""
    shell:
        "samtools faidx {input}"

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
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """

rule real_freqk_var_dedup:
    input:
        "freqk_indices/{species}.txt"
    output:
        "freqk_var_dedup/{species}.txt"
    benchmark:
        "benchmarks/freqk_var_dedup/{species}.bench"
    shell:
        """
        ./scripts/freqk var-dedup --index {input} --output {output}
        """

rule real_freqk_ref_dedup:
    input:
        index="freqk_var_dedup/{species}.txt",
        vcf=config["real_vcf"],
        fasta=config["real_fasta"]
    output:
        "freqk_ref_dedup/{species}.txt"
    benchmark:
        "benchmarks/freqk_var_dedup/{species}.bench"
    shell:
        "./scripts/freqk ref-dedup --index {input.index} --output {output} --fasta {input.fasta} --vcf {input.vcf} -k {params.k} "

rule real_freqk_count:
    input:
        index="freqk_ref_dedup/{species}.txt",
        pread1="real_fastp_results/trimmed_paired_R1_{species}.fastq.gz",
        pread2="real_fastp_results/trimmed_paired_R2_{species}.fastq.gz",
        uread1="real_fastp_results/trimmed_unpaired_R1_{species}.fastq.gz",
        uread2="real_fastp_results/trimmed_unpaired_R2_{species}.fastq.gz",
    output:
        all=temp("all_{species}.fastq"),
        counts="real_freqk_results/{species}_counts.txt",
        freqs="real_freqk_results/{species}_freqs.txt"
    benchmark:
        "benchmarks/freqk_count/{species}.bench"
    shell:
        """
        cat {input.pread1} {input.pread2} {input.uread1} {input.uread2} > {output.all}
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {output.all} --freq-output {output.freqs} --count-output {output.counts}
        """

rule real_vg_autoindex:
    input:
        fasta="ancestral_genome_results/{SID}.fasta",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz"
    output:
        dist = temp("{SID}_{PID}.dist"),
        gbz = temp("{SID}_{PID}.giraffe.gbz"),
        min = temp("{SID}_{PID}.shortread.withzip.min"),
        zip = temp("{SID}_{PID}.shortread.zipcodes")
    benchmark:
        "benchmarks/vg_autoindex/{SID}_{PID}.bench"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe -r {input.fasta} -v {input.vcf} -p {wildcards.SID}_{wildcards.PID}"

rule real_vg_giraffe_paired:
    input:
        dist = "{SID}_{PID}.dist",
        gbz = "{SID}_{PID}.giraffe.gbz",
        min = "{SID}_{PID}.shortread.withzip.min",
        zip = "{SID}_{PID}.shortread.zipcodes",
        pread1 = "fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2 = "fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/paired_{SID}_{PID}.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{SID}_{PID}.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -p -t {threads} -f {input.pread1} -f {input.pread2} > {output}
        """

rule real_vg_giraffe_unpaired:
    input:
        dist = "{SID}_{PID}.dist",
        gbz = "{SID}_{PID}.giraffe.gbz",
        min = "{SID}_{PID}.shortread.withzip.min",
        zip = "{SID}_{PID}.shortread.zipcodes",
        uread = "fastp_results/trimmed_unpaired_{read}_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/{read}_{SID}_{PID}.gam")
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{read}_{SID}_{PID}.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -p -t {threads} -f {input.uread} > {output}
        """

rule real_vg_surject:
    input:
        gbz = "{SID}_{PID}.giraffe.gbz",
        gam = "vg_giraffe_results/{read}_{SID}_{PID}.gam",
    output:
        temp("vg_surject_results/{read}_{SID}_{PID}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_surject/{read}_{SID}_{PID}.bam"
    shell:
        """
        vg surject -x {input.gbz} --progress -t {threads} -b {input.gam} > {output}
        """

rule real_samtools_sort:
    input:
        "vg_surject_results/{read}_{SID}_{PID}.bam"
    output:
        temp("vg_surject_results/sorted_{read}_{SID}_{PID}.bam")
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_sort/{read}_{SID}_{PID}.bam"
    shell:
        """
        samtools sort {input} -o {output}
        """

rule real_samtools_merge:
    input:
        expand("vg_surject_results/sorted_{read}_{{SID}}_{{PID}}.bam", read = ["paired", "R1", "R2"])
    output:
        temp("vg_surject_results/merged_{SID}_{PID}.bam")
    benchmark:
        "benchmarks/samtools_merge/{SID}_{PID}.bam"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools merge -o {output} {input}
        """

rule real_freebayes_vg:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        index="seqkit_results/ref_{SID}.fasta.fai",
        trimbam="vg_surject_results/merged_{SID}_{PID}.bam",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        tbi="slim_results/samples_{SID}_{PID}.vcf.gz.tbi"
    output:
        "freebayes_vg_results/{SID}_{PID}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/freebayes_vg/{SID}_{PID}.bench"
    params:
        n=get_pool
    shell:
        """
        freebayes -f {input.reffasta} -p {params.n} --use-best-n-alleles 2 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} 1> {output}
        """

