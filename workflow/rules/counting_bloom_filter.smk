rule counting_bloom_filter:
    input:
        "kmc_results/kmer_counts_{SID}_{PID}.txt"
    output:
        temp("cbf_results/{SID}_{PID}.txt")
    params:
        array_size = config["array_size"],
        num_hash = config["num_hash"]
    conda:
        "../envs/cbf.yaml"
    shell:
        """
        if [ ! -d "cbf_results" ]; then
            mkdir cbf_results
        fi
        
        python scripts/counting_bloom_filter.py --input {input} --output {output} --array-size {params.array_size} --num-hash {params.num_hash}
        """

