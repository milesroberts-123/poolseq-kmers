rule multiqc:
    input:
        expand("fastp_results/{SID}_0_R1R2.json", SID=parameters.index.get_level_values("ID")),
        #expand("fastp_results/{SID}_{PID}_R1R2.json", SID=two_pop_params.index.get_level_values("ID"), PID=[0,1]),
        #expand("fastp_results/{SID}_{PID}_U1.json", SID=two_pop_params.index.get_level_values("ID"), PID=[0,1]),
        #expand("fastp_results/{SID}_{PID}_U2.json", SID=two_pop_params.index.get_level_values("ID"), PID=[0,1]),
    output:
        "multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc fastp_results/
        """
