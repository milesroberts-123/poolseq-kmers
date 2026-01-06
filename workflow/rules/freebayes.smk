rule freebayes_de_novo:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        index="seqkit_results/ref_{SID}.fasta.fai",
        trimbam="bwa_results/{SID}_{PID}.bam",
    output:
        "freebayes_de_novo_results/{SID}_{PID}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/freebayes_de_novo/{SID}_{PID}.bench"
    params:
        n=get_pool
    shell:
        """
        freebayes -f {input.reffasta} -p {params.n} --use-best-n-alleles 2 --pooled-discrete {input.trimbam} 1> {output}
        """

rule freebayes_a_priori:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        index="seqkit_results/ref_{SID}.fasta.fai",
        trimbam="bwa_results/{SID}_{PID}.bam",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        tbi="slim_results/samples_{SID}_{PID}.vcf.gz.tbi"
    output:
        "freebayes_a_priori_results/{SID}_{PID}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/freebayes_a_priori/{SID}_{PID}.bench"
    params:
        n=get_pool
    shell:
        """
        freebayes -f {input.reffasta} -p {params.n} --use-best-n-alleles 2 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} 1> {output}
        """
