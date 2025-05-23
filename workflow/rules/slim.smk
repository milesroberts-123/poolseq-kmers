rule slim:
    input:
        "ancestral_genome_results/{ID}.fasta"
    output:
        temp("slim_results/{ID}.vcf"),
        temp("slim_results/{ID}.fasta")
    log:
        "logs/slim/{ID}.log"
    params:
        simtype = lookup(query="ID == '{ID}'", within=parameters, cols="simtype"),
        sigma=lookup(query="ID == '{ID}'", within=parameters, cols="sigma"),
        N=lookup(query="ID == '{ID}'", within=parameters, cols="N"),
        N1=lookup(query="ID == '{ID}'", within=parameters, cols="N1"),
        N2=lookup(query="ID == '{ID}'", within=parameters, cols="N2"),
        mg1=lookup(query="ID == '{ID}'", within=parameters, cols="mg1"),
        mg2=lookup(query="ID == '{ID}'", within=parameters, cols="mg2"),
        n=lookup(query="ID == '{ID}'", within=parameters, cols="n"),
        h=lookup(query="ID == '{ID}'", within=parameters, cols="h"),
        s=lookup(query="ID == '{ID}'", within=parameters, cols="s"),
        mu=lookup(query="ID == '{ID}'", within=parameters, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=parameters, cols="R"),
        tau=lookup(query="ID == '{ID}'", within=parameters, cols="tau"),
        qtl_mean=lookup(query="ID == '{ID}'", within=parameters, cols="qtl_mean"),
        qtl_sigma=lookup(query="ID == '{ID}'", within=parameters, cols="qtl_sigma"),
        qtl_prop=lookup(query="ID == '{ID}'", within=parameters, cols="qtl_prop"),
        optimum_mean=lookup(query="ID == '{ID}'", within=parameters, cols="optimum_mean"),
        optimum_sigma=lookup(query="ID == '{ID}'", within=parameters, cols="optimum_sigma"),
        phenotype_cutoff=lookup(query="ID == '{ID}'", within=parameters, cols="phenotype_cutoff")
    conda:
        "../envs/slim.yaml"
    shell:
        """
        #if [ -d /.singularity.d ]; then
        #    echo Singularity dectected! Activing conda env in container...
        #    mamba activate /conda-envs/5d8d7dec5540d725ad3e6cb69c16b0a2
        #fi

        if [ "{params.simtype}" == "onepop" ]; then
            slim -d ID={wildcards.ID} -d sigma={params.sigma} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/neutral.slim &> {log}
        fi

        if [ "{params.simtype}" == "twopop" ]; then
            slim -d ID={wildcards.ID} -d sigma={params.sigma} -d N1={params.N1} -d N2={params.N2} -d mg1={params.mg1} -d mg2={params.mg2} -d mu={params.mu} -d R={params.R} -d n={params.n} -d tau={params.tau} scripts/two_pop.slim &> {log}
        fi

        if [ "{params.simtype}" == "bsa" ]; then
            slim -d ID={wildcards.ID} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} -d qtl_mean={params.qtl_mean} -d qtl_sigma={params.qtl_sigma} -d qtl_prop={params.qtl_prop} -d optimum_mean={params.optimum_mean} -d optimum_sigma={params.optimum_sigma} -d phenotype_cutoff={params.phenotype_cutoff} scripts/bsa.slim &> {log}
        fi

        if [ "{params.simtype}" == "sweep" ]; then
            slim -d ID={wildcards.ID} -d h={params.h} -d s={params.s} -d sigma={params.sigma} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/sweep.slim &> {log}
        fi
        """
