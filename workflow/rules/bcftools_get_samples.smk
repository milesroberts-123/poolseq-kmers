rule bcftools_get_samples:
    input:
        samplevcf="slim_results/samples_{SID}.vcf.gz",
    output:
        popvcf=temp("slim_results/samples_{SID}_{PID}.vcf.gz"),
        tbi=temp("slim_results/samples_{SID}_{PID}.vcf.gz.tbi")
    conda:
        "../envs/bcftools.yaml"
    params:
        sammies = get_samples
    shell:
        """
        bcftools view --samples {params.sammies} -Oz -o {output.popvcf} {input}

        tabix {output.popvcf}
        """
