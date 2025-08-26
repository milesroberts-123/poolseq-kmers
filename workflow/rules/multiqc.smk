rule multiqc:
    input:
        expand("fastp_results/{SID}_0.json", SID=one_pop_params.index.get_level_values("ID")),
        #expand("fastp_results/{SID}_0_U2.json", SID=one_pop_params.index.get_level_values("ID")),
        #expand("fastp_results/{SID}_0_U1.json", SID=one_pop_params.index.get_level_values("ID")),
        #expand("fastp_results/{SID}_{PID}_R1R2.json", SID=two_pop_params.index.get_level_values("ID"), PID=[0,1]),
        #expand("fastp_results/{SID}_{PID}_U1.json", SID=two_pop_params.index.get_level_values("ID"), PID=[0,1]),
        #expand("fastp_results/{SID}_{PID}_U2.json", SID=two_pop_params.index.get_level_values("ID"), PID=[0,1]),
    output:
        "multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc.log",
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        multiqc fastp_results/
        """
