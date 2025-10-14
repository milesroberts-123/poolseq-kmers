rule vg_autoindex:
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

rule vg_giraffe_paired:
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

rule vg_giraffe_unpaired:
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

rule vg_surject:
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

rule samtools_sort:
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

rule samtools_merge:
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

rule freebayes_vg:
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
