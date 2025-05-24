rule fastp_one_pop:
    input:
        read1="iss_results/reads_{ID}_R1.fastq",
        read2="iss_results/reads_{ID}_R2.fastq",
    output:
        dpread1=temp("fastp_results/dedup_paired_R1_{ID}.fastq"),
        dpread2=temp("fastp_results/dedup_paired_R2_{ID}.fastq"),
        duread1=temp("fastp_results/dedup_unpaired_R1_{ID}.fastq"),
        duread2=temp("fastp_results/dedup_unpaired_R2_{ID}.fastq"),
        pread1=temp("fastp_results/trimmed_paired_R1_{ID}.fastq"),
        pread2=temp("fastp_results/trimmed_paired_R2_{ID}.fastq"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{ID}.fastq"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{ID}.fastq"),
        jsonR1R2="fastp_results/{ID}_R1R2.json",
        jsonU1="fastp_results/{ID}_U1.json",
        jsonU2="fastp_results/{ID}_U2.json",
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{ID}.log",
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        fastp -u {params.unqualLimit} -q {params.qualThresh} --dedup --correction -i {input.read1} -I {input.read2} -o {output.dpread1} -O {output.dpread2} --unpaired1 {output.duread1} --unpaired2 {output.duread2} &> {log}

        # trim low quality bases
        fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {output.dpread1} -I {output.dpread2} -o {output.pread1} -O {output.pread2} &> {log}
        fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU1} -i {output.duread1} -o {output.uread1} &> {log}
        fastp -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU2} -i {output.duread2} -o {output.uread2} &> {log}
        """
