rule kmc_histo:
    input:
        "../config/genomes/{species}.fna"
    output:
        histo="kmc_histo_results/{species}_{k}.histo",
        #pre=temp("counts_{species}.kmc_pre"),
        #suf=temp("counts_{species}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.species}_{wildcards.k}" ]; then
            rm -r tmp_kmc_{wildcards.species}_{wildcards.k}
        fi

        mkdir tmp_kmc_{wildcards.species}_{wildcards.k}

        # count k-mers
        # just estimate histogram only
        kmc -t{threads} -e -m9 -ci1 -cs100000 -fm -k{wildcards.k} {input} {output.histo} tmp_kmc_{wildcards.species}_{wildcards.k}

        rm -r tmp_kmc_{wildcards.species}_{wildcards.k}
        """
