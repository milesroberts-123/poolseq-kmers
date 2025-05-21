rule varscan_one_pop:
    input:
        reffasta = "seqkit_results/ref_{ID}.fasta",
        trimbam = "bwa_results/{ID}.bam"
    output:
        "varscan_results/{ID}.tsv"
    conda:
        "../envs/varscan.yaml"
    log: 
        "logs/varscan/{ID}.log"
    benchmark:
        "benchmarks/varscan/{ID}.bench"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp 1> {output} 2> {log}
        """
