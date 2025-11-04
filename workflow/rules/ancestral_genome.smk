rule ancestral_genome:
    input:
        "../config/parameters.tsv",
    output:
        "ancestral_genome_results/{ID}.fasta",
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
        Rscript scripts/random_genome.R {params.pA} {params.pC} {params.pG} {params.pT} {params.shape} {params.k} {params.L} {wildcards.ID} {params.shuffleKmers}
        """
