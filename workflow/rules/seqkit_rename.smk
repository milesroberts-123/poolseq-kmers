rule seqkit_rename:
    input:
        "unitig_caller_results/unitigs_{SID}_{PID}.fasta",
    output:
        "unitig_caller_results/unitigs_renamed_{SID}_{PID}.fasta",
    conda:
        "../envs/seqkit.yaml"
    log:
        "logs/seqkit_rename/{SID}_{PID}.log",
    shell:
        """
        seqkit rename {input} > {output}
        """
