rule smudgeplot_two_pop:
    input:
        p1="kmc_results/kmer_counts_{ID}_p1.txt",
        p2="kmc_results/kmer_counts_{ID}_p2.txt",
    output:
        "smudgeplot_results/{ID}_p1_coverages.tsv",
        "smudgeplot_results/{ID}_p1_sequences.tsv",
        "smudgeplot_results/{ID}_p2_coverages.tsv",
        "smudgeplot_results/{ID}_p2_sequences.tsv",
    conda:
        "../envs/smudgeplot.yaml"
    log:
        "logs/smudgeplot/{ID}.log",
    shell:
        """
        smudgeplot.py hetkmers -o smudgeplot_results/{wildcards.ID}_p1 --middle {input.p1} &>> {log}
        smudgeplot.py hetkmers -o smudgeplot_results/{wildcards.ID}_p2 --middle {input.p2} &>> {log}
        """
