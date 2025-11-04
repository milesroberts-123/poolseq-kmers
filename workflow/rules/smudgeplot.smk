rule smudgeplot:
    input:
        "kmc_results/kmer_counts_{SID}_{PID}.txt",
    output:
        "smudgeplot_results/{SID}_{PID}_coverages.tsv",
        "smudgeplot_results/{SID}_{PID}_sequences.tsv",
    conda:
        "../envs/smudgeplot.yaml"
    benchmark:
        "benchmarks/smudgeplot/{SID}_{PID}.bench"
    shell:
        "smudgeplot.py hetkmers -o smudgeplot_results/{wildcards.SID}_{wildcards.PID} --middle {input}"
