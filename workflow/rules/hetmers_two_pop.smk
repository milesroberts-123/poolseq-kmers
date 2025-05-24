rule hetmers_two_pop:
    input:
        p1="kmc_results/kmer_counts_{ID}_p1.txt",
        p2="kmc_results/kmer_counts_{ID}_p2.txt",
        p1p2="kmc_results/kmer_counts_{ID}_p1p2.txt",
    output:
        "hetmers_results/{ID}_p1_counts.csv",
        "hetmers_results/{ID}_p1_empirical_freqs.csv",
        "hetmers_results/{ID}_p1_bayes_states.csv",
        "hetmers_results/{ID}_p1_hashes.csv",
        "hetmers_results/{ID}_p1_seqs.csv",
        "hetmers_results/{ID}_p2_counts.csv",
        "hetmers_results/{ID}_p2_empirical_freqs.csv",
        "hetmers_results/{ID}_p2_bayes_states.csv",
        "hetmers_results/{ID}_p2_hashes.csv",
        "hetmers_results/{ID}_p2_seqs.csv",
        "hetmers_results/{ID}_p1p2_counts.csv",
        "hetmers_results/{ID}_p1p2_empirical_freqs.csv",
        "hetmers_results/{ID}_p1p2_bayes_states.csv",
        "hetmers_results/{ID}_p1p2_hashes.csv",
        "hetmers_results/{ID}_p1p2_seqs.csv",
    log:
        "logs/hetmers/{ID}.log",
    benchmark:
        "benchmarks/hetmers/{ID}.bench"
    params:
        mincount=config["mincount"],
        pool=get_pool,
        cov=lookup(query="ID == '{ID}'", within=parameters, cols="cov"),
    shell:
        """
        # create directory
        if [ ! -d "hetmers_results" ]; then
            mkdir hetmers_results
        fi

        ./scripts/hetmers --inputs {input.p1} --outputs {wildcards.ID}_p1 --coverages {params.cov} --pools {params.pool} --alphas 1 --betas 1 --minimums {params.mincount} &>> {log}

        ./scripts/hetmers --inputs {input.p2} --outputs {wildcards.ID}_p2 --coverages {params.cov} --pools {params.pool} --alphas 1 --betas 1 --minimums {params.mincount} &>> {log}

        ./scripts/hetmers --inputs {input.p1p2} --outputs {wildcards.ID}_p1p2 --coverages {params.cov} --pools {params.pool} --alphas 1 --betas 1 --minimums {params.mincount} &>> {log}

        # move output to directory
        mv {wildcards.ID}_p1_counts.csv hetmers_results/
        mv {wildcards.ID}_p1_empirical_freqs.csv hetmers_results/
        mv {wildcards.ID}_p1_bayes_states.csv hetmers_results/
        mv {wildcards.ID}_p1_hashes.csv hetmers_results/
        mv {wildcards.ID}_p1_seqs.csv hetmers_results/

        mv {wildcards.ID}_p2_counts.csv hetmers_results/
        mv {wildcards.ID}_p2_empirical_freqs.csv hetmers_results/
        mv {wildcards.ID}_p2_bayes_states.csv hetmers_results/
        mv {wildcards.ID}_p2_hashes.csv hetmers_results/
        mv {wildcards.ID}_p2_seqs.csv hetmers_results/

        mv {wildcards.ID}_p1p2_counts.csv hetmers_results/
        mv {wildcards.ID}_p1p2_empirical_freqs.csv hetmers_results/
        mv {wildcards.ID}_p1p2_bayes_states.csv hetmers_results/
        mv {wildcards.ID}_p1p2_hashes.csv hetmers_results/
        mv {wildcards.ID}_p1p2_seqs.csv hetmers_results/
        """
