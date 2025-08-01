rule bcftools_remove_ref:
    input:
        compvcf="slim_results/{ID}.vcf.gz",
        vcfidx="slim_results/{ID}.vcf.gz.tbi",
    output:
        samplevcf=temp("slim_results/samples_{ID}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/bcftools_remove_ref/{ID}.log",
    shell:
        """
        # remove reference
        # bcftools view --samples-file ^../config/ref.txt -Oz -o {output.samplevcf} {input.compvcf} &>> {log}
        bcftools view --samples ^i0 -Oz -o {output.samplevcf} {input.compvcf} &>> {log}
        """
