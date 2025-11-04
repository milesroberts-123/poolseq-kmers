rule bcftools_poolsnp:
    input:
        #ref = "seqkit_results/ref_{ID}.fasta",
        vcf="{SID}_{PID}_poolsnp_output.vcf.gz",
    output:
        tbi=temp("{SID}_{PID}_poolsnp_output.vcf.gz.tbi"),
        final="poolsnp_results/{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # unpack gzip
        gunzip {input.vcf}

        # recompress with bgzip
        bgzip $(basename {input.vcf} .gz)

        # index vcf
        tabix {input.vcf}

        # output allele depths
        bcftools view -m2 -M2 -v snps {input.vcf} | bcftools query -f '%CHROM %POS %REF %ALT [ %AD] [ %DP]\n' | sed 's:,:\t:g' > {output.final}     
        """
