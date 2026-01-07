rule bcftools_freebayes_de_novo:
    input:
        vcf="freebayes_de_novo_results/{SID}_{PID}.vcf",
    output:
        vcfgz=temp("freebayes_de_novo_results/{SID}_{PID}.vcf.gz"),
        tbi=temp("freebayes_de_novo_results/{SID}_{PID}.vcf.gz.tbi"),
        final="freebayes_de_novo_results/{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # recompress with bgzip
        bgzip {input.vcf}

        # index vcf
        tabix {output.vcfgz}

        # output allele depths
        bcftools view -m2 -M2 -v snps {output.vcfgz} | bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' > {output.final} 
        """


rule bcftools_freebayes_a_priori:
    input:
        vcf="freebayes_a_priori_results/{SID}_{PID}.vcf",
    output:
        vcfgz=temp("freebayes_a_priori_results/{SID}_{PID}.vcf.gz"),
        tbi=temp("freebayes_a_priori_results/{SID}_{PID}.vcf.gz.tbi"),
        final="freebayes_a_priori_results/{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # recompress with bgzip
        bgzip {input.vcf}

        # index vcf
        tabix {output.vcfgz}

        # output allele depths
        bcftools view -m2 -M2 -v snps {output.vcfgz} | bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' > {output.final}
        """


rule bcftools_freebayes_vg:
    input:
        vcf="freebayes_vg_results/{SID}_{PID}.vcf",
    output:
        vcfgz=temp("freebayes_vg_results/{SID}_{PID}.vcf.gz"),
        tbi=temp("freebayes_vg_results/{SID}_{PID}.vcf.gz.tbi"),
        final="freebayes_vg_results/{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # recompress with bgzip
        bgzip {input.vcf}

        # index vcf
        tabix {output.vcfgz}

        # output allele depths
        bcftools view -m2 -M2 -v snps {output.vcfgz} | bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' > {output.final}
        """
