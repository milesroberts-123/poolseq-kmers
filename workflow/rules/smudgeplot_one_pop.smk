rule smudgeplot_one_pop:
    input:
        "kmc_results/kmer_counts_{ID}.txt",
    output:
        "smudgeplot_results/{ID}_coverages.tsv",
        "smudgeplot_results/{ID}_sequences.tsv",
    conda:
        "../envs/smudgeplot.yaml"
    log:
        "logs/smudgeplot/{ID}.log",
    benchmark:
        "benchmarks/smudgeplot/{ID}.bench"
    shell:
        "smudgeplot.py hetkmers -o smudgeplot_results/{wildcards.ID} --middle {input} &> {log}"
