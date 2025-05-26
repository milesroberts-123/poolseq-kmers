rule bcftools_freebayes:
    input:
        vcf="freebayes_results/{SID}_{PID}.vcf",
    output:
        vcfgz=temp("freebayes_resuls/{SID}_{PID}.vcf.gz"),
        tbi=temp("freebayes_results/{SID}_{PID}.vcf.gz.tbi"),
        final="freebayes_results/{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/bcftools_freebayes/{SID}_{PID}.log",
    shell:
        """
        # recompress with bgzip
        bgzip {input.vcf}

        # index vcf
        tabix {input.vcf}

        # output allele depths
        bcftools view -m2 -M2 -v snps {output.vcfgz} | bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' | sed 's:,:\t:g' > {output.final}     
        """
