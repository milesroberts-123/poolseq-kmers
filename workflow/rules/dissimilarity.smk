rule dissimilarity:
    input:
        p1="kmc_results/kmer_counts_{SID}_0.txt",
        p2="kmc_results/kmer_counts_{SID}_1.txt",
    output:
        "dissimilarity_results/{SID}.txt",
    log:
        "logs/dissimilarity/{SID}.log",
    conda:
        "../envs/R.yaml"
    shell:
        """
        Rscript scripts/dissimilarity.R {input.p1} {input.p2} {output} {threads} &> {log}
        """
