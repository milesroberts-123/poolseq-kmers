rule slim:
    input:
        "ancestral_genome_results/{ID}.fasta",
    output:
        temp("slim_results/{ID}.vcf"),
        temp("slim_results/{ID}.fasta"),
    log:
        "logs/slim/{ID}.log",
    params:
        simulation=parameters.instance
        #simtype=lookup(query="ID == '{ID}'", within=parameters, cols="simtype"),
        #sigma=lookup(query="ID == '{ID}'", within=parameters, cols="sigma"),
        #N=lookup(query="ID == '{ID}'", within=parameters, cols="N"),
        #N1=lookup(query="ID == '{ID}'", within=parameters, cols="N1"),
        #N2=lookup(query="ID == '{ID}'", within=parameters, cols="N2"),
        #mg1=lookup(query="ID == '{ID}'", within=parameters, cols="mg1"),
        #mg2=lookup(query="ID == '{ID}'", within=parameters, cols="mg2"),
        #n=lookup(query="ID == '{ID}'", within=parameters, cols="n"),
        #h=lookup(query="ID == '{ID}'", within=parameters, cols="h"),
        #s=lookup(query="ID == '{ID}'", within=parameters, cols="s"),
        #mu=lookup(query="ID == '{ID}'", within=parameters, cols="mu"),
        #R=lookup(query="ID == '{ID}'", within=parameters, cols="R"),
        #tau=lookup(query="ID == '{ID}'", within=parameters, cols="tau"),
        #qtl_mean=lookup(query="ID == '{ID}'", within=parameters, cols="qtl_mean"),
        #qtl_sigma=lookup(query="ID == '{ID}'", within=parameters, cols="qtl_sigma"),
        #qtl_prop=lookup(query="ID == '{ID}'", within=parameters, cols="qtl_prop"),
        #optimum_mean=lookup(
        #    query="ID == '{ID}'", within=parameters, cols="optimum_mean"
        #),
        #optimum_sigma=lookup(
        #    query="ID == '{ID}'", within=parameters, cols="optimum_sigma"
        #),
        #phenotype_cutoff=lookup(
        #    query="ID == '{ID}'", within=parameters, cols="phenotype_cutoff"
        #),
    conda:
        "../envs/slim.yaml"
    shell:
        """
        #if [ -d /.singularity.d ]; then
        #    echo Singularity dectected! Activing conda env in container...
        #    mamba activate /conda-envs/5d8d7dec5540d725ad3e6cb69c16b0a2
        #fi

        if [ simulation["simtype"] == "onepop" ]; then
            slim -d ID={wildcards.ID} -d sigma=simulation["sigma"] -d N=simulation["N"] -d mu=simulation["mu"] -d R=simulation["R"] -d n=simulation["n"] scripts/neutral.slim &> {log}
        fi

        if [ "simulation["simtype"]" == "twopop" ]; then
            slim -d ID={wildcards.ID} -d sigma=simulation["sigma"] -d N1=simulation["N1"] -d N2=simulation["N2"] -d mg1=simulation["mg1"] -d mg2=simulation["mg2"] -d mu=simulation["mu"] -d R=simulation["R"] -d n=simulation["n"] -d tau=simulation["tau"] scripts/two_pop.slim &> {log}
        fi

        if [ "simulation["simtype"]" == "bsa" ]; then
            slim -d ID={wildcards.ID} -d N=simulation["N"] -d mu=simulation["mu"] -d R=simulation["R"] -d n=simulation["n"] -d qtl_mean=simulation["qtl_mean"] -d qtl_sigma=simulation["qtl_sigma"] -d qtl_prop=simulation["qtl_prop"] -d optimum_mean=simulation["optimum_mean"] -d optimum_sigma=simulation["optimum_sigma"] -d phenotype_cutoff=simulation["phenotype_cutoff"] scripts/bsa.slim &> {log}
        fi

        if [ "simulation["simtype"]" == "sweep" ]; then
            slim -d ID={wildcards.ID} -d h=simulation["h"] -d s=simulation["s"] -d sigma=simulation["sigma"] -d N=simulation["N"] -d mu=simulation["mu"] -d R=simulation["R"] -d n=simulation["n"] scripts/sweep.slim &> {log}
        fi
        """
