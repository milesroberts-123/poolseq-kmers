rule ancestral_genome:
    input:
        "../config/parameters.tsv",
    output:
        "ancestral_genome_results/{ID}.fasta",
    log:
        "logs/ancestral_genome/{ID}.log",
    params:
        pA=config["pA"],
        pC=config["pC"],
        pG=config["pG"],
        pT=config["pT"],
        shape=lookup(query="ID == '{ID}'", within=parameters, cols="shape"),
        L=config["L"],
        k=config["k"],
        shuffleKmers=lookup(query="ID == '{ID}'", within=parameters, cols="shuffle"),
    conda:
        "../envs/R.yaml"
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
