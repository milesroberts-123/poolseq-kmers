rule kmc:
    input:
        pread1="seqtk_results/paired_R1_{SID}_{PID}.fastq",
        pread2="seqtk_results/paired_R2_{SID}_{PID}.fastq",
        uread1="seqtk_results/unpaired_R1_{SID}_{PID}.fastq",
        uread2="seqtk_results/unpaired_R2_{SID}_{PID}.fastq",
    output:
        counts=temp("kmc_results/kmer_counts_{SID}_{PID}.txt"),
        list=temp("{SID}_{PID}.list"),
        pre=temp("tmp_counts_{SID}_{PID}.kmc_pre"),
        suf=temp("tmp_counts_{SID}_{PID}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    benchmark:
        "benchmarks/kmc/{SID}_{PID}.bench"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.SID}_{wildcards.PID}" ]; then
            rm -r tmp_kmc_{wildcards.SID}_{wildcards.PID}
        fi

        mkdir tmp_kmc_{wildcards.SID}_{wildcards.PID}

        # create file list
        echo {input.pread1} {input.pread2} {input.uread1} {input.uread2} | tr ' ' '\n' > {output.list}
        
        # count k-mers
        kmc -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} @{output.list} tmp_counts_{wildcards.SID}_{wildcards.PID} tmp_kmc_{wildcards.SID}_{wildcards.PID}

        # dump all k-mers to text file
        kmc_tools transform tmp_counts_{wildcards.SID}_{wildcards.PID} dump {output.counts} &>> {log}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.SID}_{wildcards.PID}
        """
