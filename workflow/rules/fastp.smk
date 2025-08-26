rule fastp:
    input:
        read1="iss_results/reads_{SID}_{PID}_R1.fastq",
        read2="iss_results/reads_{SID}_{PID}_R2.fastq",
    output:
        pread1=temp("fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq"),
        pread2=temp("fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{SID}_{PID}.fastq"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{SID}_{PID}.fastq"),
        jsonR1R2="fastp_results/{SID}_{PID}_R1R2.json",
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{SID}_{PID}.log",
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        # trim low quality bases
        fastp --thread {threads} --dont_eval_duplication -u {params.unqualLimit} -q {params.qualThresh} --correction -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {input.read1} -I {input.read2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2} &>> {log}
        """
