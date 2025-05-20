rule poolsnp_two_pop:
    input:
        reffasta = "seqkit_results/ref_{ID}_p1.fasta",
        trimbam_p1 = "bwa_results/{ID}_p1.bam",
        trimbam_p2 = "bwa_results/{ID}_p2.bam"
    output:
        vcf_p1 = temp("{ID}_p1_poolsnp_output.vcf.gz"),
        cov_p1 = temp("{ID}_p1_poolsnp_output-cov-0.9999.txt"),
        bs_p1 = temp("{ID}_p1_poolsnp_output_BS.txt.gz"),
        mpileup_p1 = temp("{ID}_p1.mpileup"),
        vcf_p2 = temp("{ID}_p2_poolsnp_output.vcf.gz"),
        cov_p2 = temp("{ID}_p2_poolsnp_output-cov-0.9999.txt"),
        bs_p2 = temp("{ID}_p2_poolsnp_output_BS.txt.gz"),
        mpileup_p2 = temp("{ID}_p2.mpileup")
    params:
        wd = get_wd,
        mincount = config["mincount"]
    conda:
        "../envs/poolsnp.yaml"
    resources:
        mem_mb_per_cpu=8000,
        time=239
    log:
        "logs/poolsnp/{ID}.log"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam_p1} > {output.mpileup_p1}

        samtools mpileup -f {input.reffasta} {input.trimbam_p2} > {output.mpileup_p2}

        PoolSNP.sh mpileup={params.wd}{output.mpileup_p1} reference={params.wd}{input.reffasta} names={wildcards.ID}_p1 max-cov=0.9999 min-cov={params.mincount} min-count={params.mincount} min-freq=0.01 miss-frac=0 badsites=1 allsites=0 output={params.wd}{wildcards.ID}_p1_poolsnp_output &>{log}
    
        PoolSNP.sh mpileup={params.wd}{output.mpileup_p2} reference={params.wd}{input.reffasta} names={wildcards.ID}_p2 max-cov=0.9999 min-cov={params.mincount} min-count={params.mincount} min-freq=0.01 miss-frac=0 badsites=1 allsites=0 output={params.wd}{wildcards.ID}_p2_poolsnp_output &> {log}
        """
