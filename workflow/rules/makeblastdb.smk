rule makeblastdb:
    input:
        ref="seqkit_results/ref_{SID}.fasta",
    output:
        ndb=temp("ref_{SID}.ndb"),
        nhr=temp("ref_{SID}.nhr"),
        nin=temp("ref_{SID}.nin"),
        njs=temp("ref_{SID}.njs"),
        not_dbfile=temp("ref_{SID}.not"),
        nsq=temp("ref_{SID}.nsq"),
        ntf=temp("ref_{SID}.ntf"),
        nto=temp("ref_{SID}.nto"),
    conda:
        "../envs/blast.yaml"
    log:
        "logs/makeblastdb/{SID}.log",
    shell:
        """
        makeblastdb -in {input.ref} -title $(basename {input.ref} .fasta) -dbtype nucl -out $(basename {input.ref} .fasta) &>> {log}
        """
