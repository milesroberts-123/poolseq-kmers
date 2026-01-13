rule vg_autoindex:
    input:
        fasta="ancestral_genome_results/{SID}.fasta",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
    output:
        dist=temp("{SID}_{PID}.dist"),
        gbz=temp("{SID}_{PID}.giraffe.gbz"),
        min=temp("{SID}_{PID}.shortread.withzip.min"),
        zip=temp("{SID}_{PID}.shortread.zipcodes"),
    benchmark:
        "benchmarks/vg_autoindex/{SID}_{PID}.bench"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe -r {input.fasta} -v {input.vcf} -p {wildcards.SID}_{wildcards.PID}"


rule vg_giraffe_paired:
    input:
        dist="{SID}_{PID}.dist",
        gbz="{SID}_{PID}.giraffe.gbz",
        min="{SID}_{PID}.shortread.withzip.min",
        zip="{SID}_{PID}.shortread.zipcodes",
        pread1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/paired_{SID}_{PID}.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{SID}_{PID}.bench"
    shell:
        "vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -t {threads} -f {input.pread1} -f {input.pread2} > {output}"


rule vg_giraffe_unpaired:
    input:
        dist="{SID}_{PID}.dist",
        gbz="{SID}_{PID}.giraffe.gbz",
        min="{SID}_{PID}.shortread.withzip.min",
        zip="{SID}_{PID}.shortread.zipcodes",
        uread="fastp_results/trimmed_unpaired_{read}_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/unpaired_{read}_{SID}_{PID}.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/unpaired_{read}_{SID}_{PID}.bench"
    shell:
        "vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -t {threads} -f {input.uread} > {output}"


rule vg_surject_unpaired:
    input:
        gbz="{SID}_{PID}.giraffe.gbz",
        gam="vg_giraffe_results/unpaired_{read}_{SID}_{PID}.gam",
    output:
        temp("vg_surject_results/unpaired_{read}_{SID}_{PID}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_surject/unpaired_{read}_{SID}_{PID}.bench"
    shell:
        "vg surject -x {input.gbz} -t {threads} -b {input.gam} > {output}"


rule vg_surject_paired:
    input:
        gbz="{SID}_{PID}.giraffe.gbz",
        gam="vg_giraffe_results/paired_{SID}_{PID}.gam",
    output:
        temp("vg_surject_results/paired_{SID}_{PID}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_surject/paired_{SID}_{PID}.bench"
    shell:
        "vg surject -x {input.gbz} -t {threads} -b {input.gam} > {output}"


rule samtools_sort_unpaired:
    input:
        "vg_surject_results/unpaired_{read}_{SID}_{PID}.bam",
    output:
        temp("samtools_sort_results/unpaired_{read}_{SID}_{PID}.bam"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_sort/unpaired_{read}_{SID}_{PID}.bench"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_sort_paired:
    input:
        "vg_surject_results/paired_{SID}_{PID}.bam",
    output:
        temp("samtools_sort_results/paired_{SID}_{PID}.bam"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_sort/paired_{SID}_{PID}.bench"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_merge:
    input:
        "samtools_sort_results/paired_{SID}_{PID}.bam",
        expand(
            "samtools_sort_results/unpaired_{read}_{{SID}}_{{PID}}.bam",
            read=["R1", "R2"],
        ),
    output:
        temp("samtools_merge_results/{SID}_{PID}.bam"),
    benchmark:
        "benchmarks/samtools_merge/{SID}_{PID}.bench"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "samtools merge -o {output} {input}"

rule samtools_markdup:
    input:
        "samtools_merge_results/{SID}_{PID}.bam"
    output:
        temp("samtools_markdup_results/{SID}_{PID}.bam")
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_markdup/{SID}_{PID}.bench"
    shell:
        "samtools collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} - {output}"

rule freebayes_vg:
    input:
        reffasta="ancestral_genome_results/{SID}.fasta",
        index="ancestral_genome_results/{SID}.fasta.fai",
        trimbam="samtools_markdup_results/{SID}_{PID}.bam",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        tbi="slim_results/samples_{SID}_{PID}.vcf.gz.tbi",
    output:
        "freebayes_vg_results/{SID}_{PID}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/freebayes_vg/{SID}_{PID}.bench"
    params:
        n=get_pool,
    shell:
        "freebayes -f {input.reffasta} -p {params.n} --use-best-n-alleles 2 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} 1> {output}"

rule varscan_vg:
    input:
        reffasta="ancestral_genome_results/{SID}.fasta",
        index="ancestral_genome_results/{SID}.fasta.fai",
        trimbam="samtools_markdup_results/{SID}_{PID}.bam",
    output:
        "varscan_results/{SID}_{PID}.tsv",
    conda:
        "../envs/varscan.yaml"
    benchmark:
        "benchmarks/varscan/{SID}_{PID}.bench"
    params:
        mincount=config["mincount"],
        minfreq=config["minfreq"]
    shell:
        "samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp --min-coverage {params.mincount} --min-var-freq {params.minfreq} 1> {output}"

rule poolsnp_vg:
    input:
        reffasta="ancestral_genome_results/{SID}.fasta",
        index="ancestral_genome_results/{SID}.fasta.fai",
        trimbam="samtools_markdup_results/{SID}_{PID}.bam",
    output:
        vcf=temp("{SID}_{PID}_poolsnp_output.vcf.gz"),
        cov=temp("{SID}_{PID}_poolsnp_output-cov-0.9999.txt"),
        bs=temp("{SID}_{PID}_poolsnp_output_BS.txt.gz"),
        mpileup=temp("{SID}_{PID}.mpileup"),
    params:
        wd=get_wd,
        mincount=config["mincount"],
        minfreq=config["minfreq"]
    conda:
        "../envs/poolsnp.yaml"
    benchmark:
        "benchmarks/poolsnp/{SID}_{PID}.bench"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} > {output.mpileup}

        PoolSNP.sh   \
        mpileup={params.wd}{output.mpileup} \
        reference={params.wd}{input.reffasta} \
        names={wildcards.SID} \
        max-cov=0.9999 \
        min-cov={params.mincount} \
        min-count=1 \
        min-freq={params.minfreq} \
        miss-frac=0 \
        badsites=1 \
        allsites=0 \
        output={params.wd}{wildcards.SID}_{wildcards.PID}_poolsnp_output
        """
