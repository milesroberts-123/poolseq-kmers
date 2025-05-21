rule seqkit_rename:
    input:
        "unitig_caller_results/unitigs_{ID}.fasta"
    output:
        "unitig_caller_results/unitigs_renamed_{ID}.fasta"
    conda:
        "../envs/seqkit.yaml"
    log:
        "logs/seqkit_rename/{ID}.log"   
    shell:
        """
        seqkit rename {input} > {output}
        """
