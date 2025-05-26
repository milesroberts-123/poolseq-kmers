rule freebayes:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        index="seqkit_results/ref_{SID}.fasta.fai",
        trimbam="bwa_results/{SID}_{PID}.bam",
    output:
        "freebayes_results/{SID}_{PID}.vcf",
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/freebayes/{SID}_{PID}.log",
    benchmark:
        "benchmarks/freebayes/{SID}_{PID}.bench"
    params:
        n=get_pool
    shell:
        """
        freebayes -f {input.reffasta} -p {params.n} --use-best-n-alleles 2 --pooled-discrete {input.trimbam} 1> {output} 2> {log}
        """
