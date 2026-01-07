rule multiqc:
    input:
        expand(
            "fastp_results/{SID}_0_R1R2.json",
            SID=parameters.index.get_level_values("ID"),
        ),
    output:
        "multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc fastp_results/
        """
