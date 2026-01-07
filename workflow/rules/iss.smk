def calc_nreads(wildcards):
    L = parameters.loc[parameters["ID"] == wildcards.SID, "L"]
    cov = parameters.loc[parameters["ID"] == wildcards.SID, "cov"]
    sequencer = parameters.loc[parameters["ID"] == wildcards.SID, "sequencer"]

    L = int(L.iloc[0])
    cov = int(cov.iloc[0])
    sequencer = sequencer.iloc[0]

    if sequencer == "miseq" or sequencer == "nextseq":
        nreads=(L*cov)/300
    if sequencer == "novaseq":
        nreads=(L*cov)/150
    if sequencer == "hiseq":
        nreads=(L*cov)/125

    return int(nreads)

rule iss:
    input:
        "seqkit_results/samples_{SID}_{PID}.fasta",
    output:
        temp("iss_results/reads_{SID}_{PID}_R1.fastq"),
        temp("iss_results/reads_{SID}_{PID}_R2.fastq"),
    conda:
        "../envs/iss.yaml"
    params:
        sequencer=lookup(query="ID == '{SID}'", within=parameters, cols="sequencer"),
        issseed=lookup(query="ID == '{SID}'", within=parameters, cols="issseed"),
        nreads=calc_nreads
    shell:
        """
        # read length = 300 bp
        #if [ "{params.sequencer}" == "miseq" ] || [ "{params.sequencer}" == "nextseq" ]; then
            # calculate number of reads for desired coverage level
        #    nreads=$(({params.L}*{params.cov}/300))
        #fi

        # read length 150 bp
        #if [ "{params.sequencer}" == "novaseq" ]; then
            # calculate number of reads for desired coverage level
        #    nreads=$(({params.L}*{params.cov}/150))
        #fi
        
        # read length 125 bp
        #if [ "{params.sequencer}" == "hiseq" ]; then
            # calculate number of reads for desired coverage level
        #    nreads=$(({params.L}*{params.cov}/125))
        #fi

        echo Number of reads to simulate: 
        echo {params.nreads}

        # simulate reads
        iss generate -g {input} --seed {params.issseed} --cpus {threads} --model {params.sequencer} -n {params.nreads} --abundance uniform --output iss_results/reads_{wildcards.SID}_{wildcards.PID}
        """
