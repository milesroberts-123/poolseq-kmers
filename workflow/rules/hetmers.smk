rule hetmers:
    input:
        "kmc_results/kmer_counts_{SID}_{PID}.txt",
    output:
        "hetmers_results/{SID}_{PID}_counts.csv",
        "hetmers_results/{SID}_{PID}_empirical_freqs.csv",
        "hetmers_results/{SID}_{PID}_bayes_states.csv",
        "hetmers_results/{SID}_{PID}_hashes.csv",
        "hetmers_results/{SID}_{PID}_seqs.csv",
    benchmark:
        "benchmarks/hetmers/{SID}_{PID}.bench"
    params:
        mincount=config["mincount"],
        pool=get_pool,
        cov=lookup(query="ID == '{SID}'", within=parameters, cols="cov"),
    shell:
        """
        # create directory
        if [ ! -d "hetmers_results" ]; then
            mkdir hetmers_results
        fi

        ./scripts/hetmers --inputs {input} --outputs {wildcards.SID}_{wildcards.PID} --coverages {params.cov} --pools {params.pool} --alphas 1 --betas 1 --minimums {params.mincount}

        # move output to directory
        mv {wildcards.SID}_{wildcards.PID}_counts.csv hetmers_results/
        mv {wildcards.SID}_{wildcards.PID}_empirical_freqs.csv hetmers_results/
        mv {wildcards.SID}_{wildcards.PID}_bayes_states.csv hetmers_results/
        mv {wildcards.SID}_{wildcards.PID}_hashes.csv hetmers_results/
        mv {wildcards.SID}_{wildcards.PID}_seqs.csv hetmers_results/
        """
