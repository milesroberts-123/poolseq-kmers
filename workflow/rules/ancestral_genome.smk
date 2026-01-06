rule ancestral_genome:
    input:
        "../config/parameters.tsv",
    output:
        "ancestral_genome_results/{ID}.fasta",
    params:
        pA=lookup(query="ID == '{ID}'", within=parameters, cols="pA"),
        pC=lookup(query="ID == '{ID}'", within=parameters, cols="pC"),
        pG=lookup(query="ID == '{ID}'", within=parameters, cols="pG"),
        pT=lookup(query="ID == '{ID}'", within=parameters, cols="pT"),
        shape=lookup(query="ID == '{ID}'", within=parameters, cols="shape"),
        L=lookup(query="ID == '{ID}'", within=parameters, cols="L"),
        k=lookup(query="ID == '{ID}'", within=parameters, cols="k"),
        shuffleKmers=lookup(query="ID == '{ID}'", within=parameters, cols="shuffle"),
    conda:
        "../envs/R.yaml"
    shell:
        """
        Rscript scripts/random_genome.R {params.pA} {params.pC} {params.pG} {params.pT} {params.shape} {params.k} {params.L} {wildcards.ID} {params.shuffleKmers}
        """
