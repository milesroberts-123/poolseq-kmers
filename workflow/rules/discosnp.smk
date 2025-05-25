rule discosnp:
    input:
        pread1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
        uread1="fastp_results/trimmed_unpaired_R1_{SID}_{PID}.fastq",
        uread2="fastp_results/trimmed_unpaired_R2_{SID}_{PID}.fastq",
        ref="seqkit_results/ref_{SID}.fasta",
        amb="seqkit_results/ref_{SID}.fasta.amb",
        ann="seqkit_results/ref_{SID}.fasta.ann",
        bwt="seqkit_results/ref_{SID}.fasta.bwt",
        pac="seqkit_results/ref_{SID}.fasta.pac",
        sa="seqkit_results/ref_{SID}.fasta.sa",
    output:
        #tmpread = temp("tmp_read_set_{SID}.fastq"),
        fasta=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + "_D_100_P_3_b_0_coherent.fa"
        ),
        fof=temp("fof_{SID}_{PID}.txt"),
        fof_reads=temp("fof_reads_{SID}_{PID}.txt"),
        vcf=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + "_D_100_P_3_b_0_coherent.vcf"
        ),
        covh5=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + "_cov.h5"
        ),
        h5=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + ".h5"
        ),
        sam=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + "_D_100_P_3_b_0_coherentBWA_MEM.sam"
        ),
        igv=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + "_D_100_P_3_b_0_coherent_for_IGV.vcf"
        ),
        uncofa=temp(
            "discoRes_{SID}_{PID}_k_"
            + str(config["k"])
            + "_c_"
            + str(config["mincount"])
            + "_D_100_P_3_b_0_uncoherent.fa"
        ),
        corres=temp("discoRes_{SID}_{PID}_read_files_correspondance.txt"),
    conda:
        "../envs/discosnp.yaml"
    log:
        "logs/discosnp/{SID}_{PID}.log",
    benchmark:
        "benchmarks/discosnp/{SID}_{PID}.bench"
    params:
        prefix="discoRes_{SID}_{PID}",
        mincount=config["mincount"],
        k=config["k"],
    priority: 100
    shell:
        """
        # create temp directory
        if [ -d "tmp_discosnp_{wildcards.SID}" ]; then
            rm -r tmp_discosnp_{wildcards.SID}
        fi

        mkdir tmp_discosnp_{wildcards.SID}

        # create file of files
        echo "../{output.fof_reads}" > {output.fof}

        # check each file for being empty, use only non-empty files
        if [ -s {input.pread1} ]; then
            echo "../{input.pread1}" >> {output.fof_reads}
        fi

        if [ -s {input.pread2} ]; then
            echo "../{input.pread2}" >> {output.fof_reads}
        fi

        if [ -s {input.uread1} ]; then
            echo "../{input.uread1}" >> {output.fof_reads}
        fi

        if [ -s {input.uread2} ]; then
            echo "../{input.uread2}" >> {output.fof_reads}
        fi

        # run discosnp, with results for mapping SNPs to reference
        cd tmp_discosnp_{wildcards.SID}

        run_discoSnp++.sh -r ../{output.fof} -c {params.mincount} -k {params.k} -G ../{input.ref} -p {params.prefix} &> ../{log}

        # move output from temp directory
        mv {params.prefix}* ..

        # delete temporary directory
        cd ..
        rm -r tmp_discosnp_{wildcards.SID}
        """
