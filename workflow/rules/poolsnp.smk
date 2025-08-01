rule poolsnp:
    input:
        reffasta="seqkit_results/ref_{SID}.fasta",
        trimbam="bwa_results/{SID}_{PID}.bam",
    output:
        vcf=temp("{SID}_{PID}_poolsnp_output.vcf.gz"),
        cov=temp("{SID}_{PID}_poolsnp_output-cov-0.9999.txt"),
        bs=temp("{SID}_{PID}_poolsnp_output_BS.txt.gz"),
        mpileup=temp("{SID}_{PID}.mpileup"),
    params:
        wd=get_wd,
        mincount=config["mincount"],
    conda:
        "../envs/poolsnp.yaml"
    benchmark:
        "benchmarks/poolsnp/{SID}_{PID}.bench"
    log:
        "logs/poolsnp/{SID}_{PID}.log",
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} > {output.mpileup}

        PoolSNP.sh   \
        mpileup={params.wd}{output.mpileup} \
        reference={params.wd}{input.reffasta} \
        names={wildcards.SID} \
        max-cov=0.9999 \
        min-cov={params.mincount} \
        min-count={params.mincount} \
        min-freq=0.01 \
        miss-frac=0 \
        badsites=1 \
        allsites=0 \
        output={params.wd}{wildcards.SID}_{wildcards.PID}_poolsnp_output &> {log}
        """
