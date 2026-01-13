rule exact_kmc_histo:
    input:
        "ancestral_genome_results/{ID}.fasta",
    output:
        histo="exact_kmc_histo_results/{ID}.histo",
        pre=temp("counts_{ID}.kmc_pre"),
        suf=temp("counts_{ID}.kmc_suf"),
    conda:
        "../envs/kmc.yaml"
    params:
        k=lookup(query="ID == '{ID}'", within=parameters, cols="k"),
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}" ]; then
            rm -r tmp_kmc_{wildcards.ID}
        fi

        mkdir tmp_kmc_{wildcards.ID}

        # count k-mers
        kmc -t{threads} -m9 -ci1 -cs100000 -fm -k{params.k} {input} counts_{wildcards.ID} tmp_kmc_{wildcards.ID}

        # convert to histogram
        kmc_tools transform counts_{wildcards.ID} histogram {output.histo}

        # rm tmp dir
        rm -r tmp_kmc_{wildcards.ID}
        """
