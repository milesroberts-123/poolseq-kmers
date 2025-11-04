rule blast:
    input:
        ref="seqkit_results/ref_{SID}.fasta",
        ndb="ref_{SID}.ndb",
        nhr="ref_{SID}.nhr",
        nin="ref_{SID}.nin",
        njs="ref_{SID}.njs",
        not_dbfile="ref_{SID}.not",
        nsq="ref_{SID}.nsq",
        ntf="ref_{SID}.ntf",
        nto="ref_{SID}.nto",
        unitigs="unitig_caller_results/unitigs_renamed_{SID}_{PID}.fasta",
    output:
        alignments="blast_results/{SID}_{PID}.txt",
    conda:
        "../envs/blast.yaml"
    params:
        blastEvalue=config["blastEvalue"],
    shell:
        """
        blastn -query {input.unitigs} -db $(basename {input.ref} .fasta) -out {output.alignments} -max_target_seqs 1 -evalue {params.blastEvalue} -outfmt 6
        """
