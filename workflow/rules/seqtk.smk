rule seqtk:
    input:
        pread1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
        uread1="fastp_results/trimmed_unpaired_R1_{SID}_{PID}.fastq",
        uread2="fastp_results/trimmed_unpaired_R2_{SID}_{PID}.fastq"
    output:
        pread1="seqtk_results/paired_R1_{SID}_{PID}.fastq",
        pread2="seqtk_results/paired_R2_{SID}_{PID}.fastq",
        uread1="seqtk_results/unpaired_R1_{SID}_{PID}.fastq",
        uread2="seqtk_results/unpaired_R2_{SID}_{PID}.fastq",
    conda:
        "../envs/seqtk.yaml"
    params:
        qualThresh=config["qualThresh"],
    shell:
        """
        seqtk seq -q {params.qualThresh} -n N {input.pread1} > {output.pread1}
        seqtk seq -q {params.qualThresh} -n N {input.pread2} > {output.pread2}
        seqtk seq -q {params.qualThresh} -n N {input.uread1} > {output.uread1}
        seqtk seq -q {params.qualThresh} -n N {input.uread2} > {output.uread2}
        """

