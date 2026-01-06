rule compress:
    input:
        "slim_results/{ID}.vcf",
    output:
        temp("slim_results/{ID}.vcf.gz"),
        temp("slim_results/{ID}.vcf.gz.tbi"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bgzip {input}
        tabix {input}.gz
        """
