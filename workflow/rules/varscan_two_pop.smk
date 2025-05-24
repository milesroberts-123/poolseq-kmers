rule varscan_two_pop:
    input:
        reffasta="seqkit_results/ref_{ID}_p1.fasta",
        trimbam_p2="bwa_results/{ID}_p1.bam",
        trimbam_p1="bwa_results/{ID}_p2.bam",
    output:
        cp1="varscan_results/{ID}_p1.tsv",
        cp2="varscan_results/{ID}_p2.tsv",
    conda:
        "../envs/varscan.yaml"
    log:
        "logs/varscan/{ID}.log",
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam_p1} | varscan pileup2snp 1> {output.cp1} 2> {log}

        samtools mpileup -f {input.reffasta} {input.trimbam_p2} | varscan pileup2snp 1> {output.cp2} 2> {log}
        """
