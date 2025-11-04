rule bcftools_discosnp:
    input:
        ref="seqkit_results/ref_{SID}.fasta",
        fai="seqkit_results/ref_{SID}.fasta.fai",
        vcf="discoRes_{SID}_{PID}_k_"
        + str(config["k"])
        + "_c_"
        + str(config["mincount"])
        + "_D_100_P_3_b_0_coherent.vcf",
    output:
        header=temp("discoRes_header_{SID}_{PID}.vcf"),
        bgzip=temp("discoRes_sorted_{SID}_{PID}.vcf.gz"),
        tbi=temp("discoRes_sorted_{SID}_{PID}.vcf.gz.tbi"),
        final="disco_results/{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # reheader
        # also, fix type definition for SNP < k bp apart
        bcftools reheader --fai {input.fai} {input.vcf} | sed 's:TySNP:Ty=SNP:g' > {output.header}

        # sort vcf
        bcftools sort {output.header} > discoRes_sorted_{wildcards.SID}_{wildcards.PID}.vcf

        # index vcf
        bgzip discoRes_sorted_{wildcards.SID}_{wildcards.PID}.vcf
        tabix {output.bgzip}

        # output allele depths
        bcftools query -f '%CHROM %POS %REF %ALT [ %AD]\n' {output.bgzip} | sed 's:,:\t:g' > {output.final}
        """
