rule kmc_two_pop:
    input:
                pread1_p1 = "fastp_results/trimmed_paired_R1_{ID}_p1.fastq",
                pread2_p1 = "fastp_results/trimmed_paired_R2_{ID}_p1.fastq",
                uread1_p1 = "fastp_results/trimmed_unpaired_R1_{ID}_p1.fastq",
                uread2_p1 = "fastp_results/trimmed_unpaired_R2_{ID}_p1.fastq",
                pread1_p2 = "fastp_results/trimmed_paired_R1_{ID}_p2.fastq",
                pread2_p2 = "fastp_results/trimmed_paired_R2_{ID}_p2.fastq",
                uread1_p2 = "fastp_results/trimmed_unpaired_R1_{ID}_p2.fastq",
                uread2_p2 = "fastp_results/trimmed_unpaired_R2_{ID}_p2.fastq"
    output:
        p1=temp("kmc_results/kmer_counts_{ID}_p1.txt"),
        p2=temp("kmc_results/kmer_counts_{ID}_p2.txt"),
        tmp_R1_p1_pre=temp("tmp_R1_{ID}_p1.kmc_pre"),
        tmp_R1_p1_suf=temp("tmp_R1_{ID}_p1.kmc_suf"),
        tmp_R2_p1_pre=temp("tmp_R2_{ID}_p1.kmc_pre"),
        tmp_R2_p1_suf=temp("tmp_R2_{ID}_p1.kmc_suf"),
        tmp_u_R1_p1_pre=temp("tmp_u_R1_{ID}_p1.kmc_pre"),
        tmp_u_R1_p1_suf=temp("tmp_u_R1_{ID}_p1.kmc_suf"),
        tmp_u_R2_p1_pre=temp("tmp_u_R2_{ID}_p1.kmc_pre"),
        tmp_u_R2_p1_suf=temp("tmp_u_R2_{ID}_p1.kmc_suf"),
        union_R1_R2_p1_pre=temp("union_R1_R2_{ID}_p1.kmc_pre"),
        union_R1_R2_p1_suf=temp("union_R1_R2_{ID}_p1.kmc_suf"),
        union_R1_R2_u1_p1_pre=temp("union_R1_R2_u1_{ID}_p1.kmc_pre"),
        union_R1_R2_u1_p1_suf=temp("union_R1_R2_u1_{ID}_p1.kmc_suf"),
        union_R1_R2_u1_u2_p1_pre=temp("union_R1_R2_u1_u2_{ID}_p1.kmc_pre"),
        union_R1_R2_u1_u2_p1_suf=temp("union_R1_R2_u1_u2_{ID}_p1.kmc_suf"),
        tmp_R1_p2_pre=temp("tmp_R1_{ID}_p2.kmc_pre"),
        tmp_R1_p2_suf=temp("tmp_R1_{ID}_p2.kmc_suf"),
        tmp_R2_p2_pre=temp("tmp_R2_{ID}_p2.kmc_pre"),
        tmp_R2_p2_suf=temp("tmp_R2_{ID}_p2.kmc_suf"),
        tmp_u_R1_p2_pre=temp("tmp_u_R1_{ID}_p2.kmc_pre"),
        tmp_u_R1_p2_suf=temp("tmp_u_R1_{ID}_p2.kmc_suf"),
        tmp_u_R2_p2_pre=temp("tmp_u_R2_{ID}_p2.kmc_pre"),
        tmp_u_R2_p2_suf=temp("tmp_u_R2_{ID}_p2.kmc_suf"),
        union_R1_R2_p2_pre=temp("union_R1_R2_{ID}_p2.kmc_pre"),
        union_R1_R2_p2_suf=temp("union_R1_R2_{ID}_p2.kmc_suf"),
        union_R1_R2_u1_p2_pre=temp("union_R1_R2_u1_{ID}_p2.kmc_pre"),
        union_R1_R2_u1_p2_suf=temp("union_R1_R2_u1_{ID}_p2.kmc_suf"),
        union_R1_R2_u1_u2_p2_pre=temp("union_R1_R2_u1_u2_{ID}_p2.kmc_pre"),
        union_R1_R2_u1_u2_p2_suf=temp("union_R1_R2_u1_u2_{ID}_p2.kmc_suf"),
    conda:
        "../envs/kmc.yaml"
    log: 
        "logs/kmc/{ID}.log"
    params:
        mincount = config["mincount"],
        maxcount = config["maxcount"],
        k = config["k"]
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}_p1" ]; then
            rm -r tmp_kmc_{wildcards.ID}_p1
        fi

        if [ -d "tmp_kmc_{wildcards.ID}_p2" ]; then
            rm -r tmp_kmc_{wildcards.ID}_p2
        fi

        mkdir tmp_kmc_{wildcards.ID}_p1

        mkdir tmp_kmc_{wildcards.ID}_p2

        # count k-mers
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread1_p1} tmp_R1_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread2_p1} tmp_R2_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.uread1_p1} tmp_u_R1_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.uread2_p1} tmp_u_R2_{wildcards.ID}_p1 tmp_kmc_{wildcards.ID}_p1 &>> {log}

        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread1_p2} tmp_R1_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread2_p2} tmp_R2_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.uread1_p2} tmp_u_R1_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}
        kmc -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.uread2_p2} tmp_u_R2_{wildcards.ID}_p2 tmp_kmc_{wildcards.ID}_p2 &>> {log}

        # combine k-mer counts into one database
        kmc_tools simple tmp_R1_{wildcards.ID}_p1 tmp_R2_{wildcards.ID}_p1 union union_R1_R2_{wildcards.ID}_p1 &>> {log}
        kmc_tools simple union_R1_R2_{wildcards.ID}_p1 tmp_u_R1_{wildcards.ID}_p1 union union_R1_R2_u1_{wildcards.ID}_p1 &>> {log}
        kmc_tools simple union_R1_R2_u1_{wildcards.ID}_p1 tmp_u_R2_{wildcards.ID}_p1 union union_R1_R2_u1_u2_{wildcards.ID}_p1 &>> {log}

        kmc_tools simple tmp_R1_{wildcards.ID}_p2 tmp_R2_{wildcards.ID}_p2 union union_R1_R2_{wildcards.ID}_p2 &>> {log}
        kmc_tools simple union_R1_R2_{wildcards.ID}_p2 tmp_u_R1_{wildcards.ID}_p2 union union_R1_R2_u1_{wildcards.ID}_p2 &>> {log}
        kmc_tools simple union_R1_R2_u1_{wildcards.ID}_p2 tmp_u_R2_{wildcards.ID}_p2 union union_R1_R2_u1_u2_{wildcards.ID}_p2 &>> {log}

        # dump all k-mers to text file
        kmc_tools transform union_R1_R2_u1_u2_{wildcards.ID}_p1 dump {output.p1} &>> {log}

        kmc_tools transform union_R1_R2_u1_u2_{wildcards.ID}_p2 dump {output.p2} &>> {log}

        # delete tmp directory
        rm -r tmp_kmc_{wildcards.ID}_p1
        rm -r tmp_kmc_{wildcards.ID}_p2
        """
