rule slim:
    input:
        "ancestral_genome_results/{ID}.fasta",
    output:
        temp("slim_results/{ID}.vcf"),
        temp("slim_results/{ID}.fasta"),
    params:
        simtype=lookup(query="ID == '{ID}'", within=parameters, cols="simtype"),
        slimseed=lookup(query="ID == '{ID}'", within=parameters, cols="slimseed"),
        N=lookup(query="ID == '{ID}'", within=parameters, cols="N"),
        n=lookup(query="ID == '{ID}'", within=parameters, cols="n"),
        h=lookup(query="ID == '{ID}'", within=parameters, cols="h"),
        s=lookup(query="ID == '{ID}'", within=parameters, cols="s"),
        mu=lookup(query="ID == '{ID}'", within=parameters, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=parameters, cols="R"),
    conda:
        "../envs/slim.yaml"
    shell:
        """
        if [ "{params.simtype}" == "onepop" ]; then
            slim -d ID={wildcards.ID} -d SLIMSEED={params.slimseed} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/neutral.slim
        fi

        if [ "{params.simtype}" == "sweep" ]; then
            slim -d ID={wildcards.ID} -d SLIMSEED={params.slimseed} -d h={params.h} -d s={params.s} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/sweep.slim
        fi
        """
