rule varscan:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        trimbam="bwa_results/{SID}_{PID}.bam",
    output:
        "varscan_results/{SID}_{PID}.tsv",
    conda:
        "../envs/varscan.yaml"
    log:
        "logs/varscan/{SID}_{PID}.log",
    benchmark:
        "benchmarks/varscan/{SID}_{PID}.bench"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp 1> {output} 2> {log}
        """
