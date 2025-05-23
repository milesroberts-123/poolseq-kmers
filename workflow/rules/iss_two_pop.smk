rule iss_two_pop:
    input:
        p1="seqkit_results/samples_{ID}_p1.fasta",
        p2="seqkit_results/samples_{ID}_p2.fasta"
    output:
        temp("iss_results/reads_{ID}_p1_R1.fastq"),
        temp("iss_results/reads_{ID}_p1_R2.fastq"),
        temp("iss_results/reads_{ID}_p2_R1.fastq"),
        temp("iss_results/reads_{ID}_p2_R2.fastq")
    conda:
        "../envs/iss.yaml"
    log:
        "logs/iss/{ID}.log"
    params:
        L = config["L"],
        cov = lookup(query="ID == '{ID}'", within=parameters, cols="cov"),
        sequencer = lookup(query="ID == '{ID}'", within=parameters, cols="sequencer")
    shell:
        """
        # read length = 300 bp
        if [ "{params.sequencer}" == "miseq" ] || [ "{params.sequencer}" == "nextseq" ]; then
            # calculate number of reads for desired coverage level
            nreads=$(({params.L}*{params.cov}/300))
        fi

        # read length 150 bp
        if [ "{params.sequencer}" == "novaseq" ]; then
            # calculate number of reads for desired coverage level
            nreads=$(({params.L}*{params.cov}/150))
        fi
        
        # read length 125 bp
        if [ "{params.sequencer}" == "hiseq" ]; then
            # calculate number of reads for desired coverage level
            nreads=$(({params.L}*{params.cov}/125))
        fi

        iss generate -g {input.p1} --cpus {threads} --model {params.sequencer} -n $nreads --abundance uniform --output iss_results/reads_{wildcards.ID}_p1 &> {log}
        iss generate -g {input.p2} --cpus {threads} --model {params.sequencer} -n $nreads --abundance uniform --output iss_results/reads_{wildcards.ID}_p2 &> {log}
        """
