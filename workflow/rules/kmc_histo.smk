rule kmc_histo:
    input:
        "../config/genomes/{species}.fna"
    output:
        histo="kmc_histo_results/{species}.histo",
        pre=temp("counts_{species}.kmc_pre"),
        suf=temp("counts_{species}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.species}" ]; then
            rm -r tmp_kmc_{wildcards.species}
        fi

        mkdir tmp_kmc_{wildcards.species}

        # count k-mers
        kmc -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread1} counts_{wildcards.species} tmp_kmc_{wildcards.species} &>> {log}

        # create histogram
        kmc_tools transform counts_{wildcards.species} histogram {output.histo}
        """
