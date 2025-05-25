rule bcftools_get_samples:
    input:
        samplevcf="slim_results/samples_{SID}.vcf.gz",
    output:
        popvcf=temp("slim_results/samples_{SID}_{PID}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    params:
        sammies = get_samples
    log:
        "logs/bcftools_get_samples/{SID}_{PID}.log",
    shell:
        """
        echo {params.sammies} &>> {log}
        bcftools view --samples {params.sammies} -Oz -o {output.popvcf} {input} &>> {log}
        """
