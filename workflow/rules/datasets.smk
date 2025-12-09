checkpoint datasets:
    output:
        directory("ncbi_datasets_results/")
    params:
        taxa = config["ncbi_taxa"]
    conda:
        "../envs/datasets.yaml"
    shell:
        r"""
        if [ -d "contam_genomes" ]; then
            rm -r contam_genomes/
        fi

        # get all contaminating genomes
        datasets download genome taxon {params.taxa} --reference --dehydrated --filename contam.zip

        # unpack metadata
        unzip contam.zip -d contam_genomes

        # download genomes based on metadata
        datasets rehydrate --max-workers {threads} --directory contam_genomes/

        # search directory for all genomes and copy them into one file
        find contam_genomes/ncbi_dataset/data/ -type f -name '*.fna' -exec cp {{}} \; > {output}

        # clean up
        # rm -r contam_genomes/
        """

rule kmc_histo:
    input:
        "ncbi_datasets_results/{species}.fna"
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

