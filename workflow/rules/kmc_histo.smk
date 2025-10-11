rule kmc_histo:
    input:
        "../config/genomes/{species}.fna"
    output:
        histo="kmc_histo_results/{species}.histo",
        #pre=temp("counts_{species}.kmc_pre"),
        #suf=temp("counts_{species}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    params:
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.species}" ]; then
            rm -r tmp_kmc_{wildcards.species}
        fi

        mkdir tmp_kmc_{wildcards.species}

        # count k-mers
        # just estimate histogram only
        kmc -t{threads} -e -m9 -ci1 -cs100000 -fm -k{params.k} {input} {output.histo} tmp_kmc_{wildcards.species}

        rm -r tmp_kmc_{wildcards.species}
        """
