def sample_start(wildcards):
    # get sample size
    n = parameters.loc[parameters["ID"] == wildcards.SID, "n"]
    n = int(n.iloc[0])

    # get population
    p = int(wildcards.PID)

    return 1+n*p

def sample_end(wildcards):
    # get sample size
    n = parameters.loc[parameters["ID"] == wildcards.SID, "n"]
    n = int(n.iloc[0])

    # get population
    p = int(wildcards.PID)

    return n*(p+1)

rule seqkit_get_samples:
    input:
        slimfasta="slim_results/{SID}.fasta",
    output:
        tempsamplefasta=temp("seqkit_results/temp_{SID}_{PID}.fasta"),
        pop1="seqkit_results/samples_{SID}_{PID}.fasta",
    params:
        start = sample_start,
        end = sample_end
    conda:
        "../envs/seqkit.yaml"
    log:
        "logs/seqkit_get_samples/{SID}_{PID}.log",
    shell:
        """
        # create file with reference genome removed
        seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.tempsamplefasta}

        echo {params.start} &>> {log}
        echo {params.end} &>> {log}

        # get group of n individuals
        seqkit range -r {params.start}:{params.end} {output.tempsamplefasta} > {output.pop1}        
        """
