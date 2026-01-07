rule bcftools_fill_tags:
    input:
        "slim_results/samples_{SID}_{PID}.vcf.gz",
    output:
        filledvcf=temp("slim_results/filled_{SID}_{PID}.vcf.gz"),
        allelefreq="slim_results/allele_freqs_{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # calculate allele frequencies
        bcftools +fill-tags {input} -Oz -o {output.filledvcf}
        bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' -o {output.allelefreq} {output.filledvcf}
        """
