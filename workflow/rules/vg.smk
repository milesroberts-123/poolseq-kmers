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
        "vg autoindex -w giraffe -g {input} -p {wildcards.SID}_{wildcards.PID}"

rule vg_giraffe:
    input:
        dist = "{SID}_{PID}.dist",
        gbz = "{SID}_{PID}.giraffe.gbz",
        min = "{SID}_{PID}.shortread.withzip.min",
        zip = "{SID}_{PID}.shortread.zipcodes",
        read1 = "fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        read2 = "fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/{SID}_{PID}.gam")
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{SID}_{PID}.bench"
    shell:
        """
        vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -p -f {input.read1} -f {input.read2} > {output}
        """

rule vg_surject:
    input:
        gbz = "{SID}_{PID}.giraffe.gbz",
        gam = "vg_giraffe_results/{SID}_{PID}.gam"
    output:
        temp("vg_surject_results/{SID}_{PID}.bam")
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_surject/{SID}_{PID}.bam"
    shell:
        "vg surject -x {input.gbz} --progress -t {threads} -b {input.gam} > {ouput}"

rule samtools_sort:
    input:
        "vg_surject_results/{SID}_{PID}.bam"
    output:
        temp("vg_surject_results/sorted_{SID}_{PID}.bam")
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_sort/{SID}_{PID}.bam"
    shell:
        "samtools sort {input} -o {output}"

rule freebayes_vg:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        index="seqkit_results/ref_{SID}.fasta.fai",
        trimbam="vg_surject_results/sorted_{SID}_{PID}.bam",
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
