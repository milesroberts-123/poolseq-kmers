rule ancestral_genome:
    input:
        "../config/parameters.tsv"
    output:
        "ancestral_genome_results/{ID}.fasta",
        #"power_law_{ID}.jpg"
    log:
        "logs/ancestral_genome/{ID}.log"
    params:
        pA=get_pA,
        pC=get_pC,
        pG=get_pG,
        pT=get_pT,
        shape=get_shape,
        L=get_L,
        k=config["k"],
        shuffleKmers=get_shuffle
    conda:
        "../envs/R.yaml"
    threads: 1
    resources:
        mem_mb_per_cpu=8000,
        time=239
    shell:
        """
        #if [ -d /.singularity.d ]; then
        #    echo Singularity detected! Using conda env in container
        #    conda init
        #    source ~/.bashrc
        #    conda activate /conda-envs/3c74495c6c8f8b13c8f5e41bfa52b11d
        #fi

        # check for singularity
        #./scripts/check_for_singularity.sh /conda-envs/3c74495c6c8f8b13c8f5e41bfa52b11d &>> {log}

        Rscript scripts/random_genome.R {params.pA} {params.pC} {params.pG} {params.pT} {params.shape} {params.k} {params.L} {wildcards.ID} {params.shuffleKmers} &>> {log}
        """
