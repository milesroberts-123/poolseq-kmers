rule iss:
    input:
        "seqkit_results/samples_{SID}_{PID}.fasta",
    output:
        temp("iss_results/reads_{SID}_{PID}_R1.fastq"),
        temp("iss_results/reads_{SID}_{PID}_R2.fastq"),
    conda:
        "../envs/iss.yaml"
    params:
        L=lookup(query="ID == '{SID}'", within=parameters, cols="L"),,
        cov=lookup(query="ID == '{SID}'", within=parameters, cols="cov"),
        sequencer=lookup(query="ID == '{SID}'", within=parameters, cols="sequencer"),
        issseed=lookup(query="ID == '{SID}'", within=parameters, cols="issseed")
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

        # simulate reads
        iss generate -g {input} --seed {params.issseed} --cpus {threads} --model {params.sequencer} -n $nreads --abundance uniform --output iss_results/reads_{wildcards.SID}_{wildcards.PID}
        """
