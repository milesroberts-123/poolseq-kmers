def check_simtype(wildcards):
    simtype = parameters.loc[parameters["ID"] == wildcards.SID, "simtype"]
    simtype = str(simtype.iloc[0])
    if (simtype == "onepop") or (simtype == "sweep"):
        return "kmc_results/kmer_counts_" + wildcards.SID + "_0.txt"
    if (simtype == "twopop") or (simtype == "bsa"):
        return ["kmc_results/kmer_counts_" + wildcards.SID + "_0.txt", "kmc_results/kmer_counts_" + wildcards.SID + "_1.txt"]

rule dissimilarity:
    input:
        check_simtype
    output:
        "dissimilarity_results/{SID}.txt",
    log:
        "logs/dissimilarity/{SID}.log",
    conda:
        "../envs/R.yaml"
    params:
        simtype=lookup(query="ID == '{SID}'", within=parameters, cols="simtype")
    shell:
        """
        if [ {params.simtype} == "twopop"] || [ {params.simtype} == "bsa"]; then
            Rscript scripts/dissimilarity.R {input} {output} {threads} &> {log}
        fi

        if [ {params.simtype} == "onepop"] || [ {params.simtype} == "sweep"]; then
            touch {output}
        fi
        """
