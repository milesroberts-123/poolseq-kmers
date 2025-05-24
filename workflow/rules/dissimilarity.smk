rule dissimilarity:
    input:
        p1="kmer_counts_{ID}_p1.txt",
        p2="kmer_counts_{ID}_p2.txt",
    output:
        "dissimilarity_{ID}.txt",
    log:
        "logs/dissimilarity/{ID}.log",
    conda:
        "../envs/R.yaml"
    shell:
        """
        Rscript scripts/dissimilarity.R {input.p1} {input.p2} {output} {threads} &> {log}
                """
