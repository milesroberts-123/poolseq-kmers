rule bedtools:
    input:
        fasta = "ref_{ID}.fasta",
        genome="../config/ref.genome",
        bed = "../config/mask.bed"
    output:
        masked_ref = "ref_masked_{ID}.fasta",
        shuf_bed = temp("shuf_{ID}.bed")
    threads: 1
    resources:
        mem_mb_per_cpu=8000,
        time=239
    conda:
        "../envs/bedtools.yaml"
    log: 
        "logs/bedtools/{ID}.log"
    shell:
        """
        # randomly place masks along genome
        bedtools shuffle -noOverlapping -i {input.bed} -g {input.genome} > {output.shuf_bed}

        # mask reference genome
        bedtools maskfasta -fi {input.fasta} -bed {output.shuf_bed} -fo {output.masked_ref}
        """
