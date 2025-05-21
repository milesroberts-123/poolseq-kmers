rule compress:
    input:
        "slim_results/{ID}.vcf",
    output:
        temp("slim_results/{ID}.vcf.gz"),
        temp("slim_results/{ID}.vcf.gz.tbi")
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/compress/{ID}.log"
    shell:
        """
        bgzip {input} &> {log}
        tabix {input}.gz &> {log}
        """
