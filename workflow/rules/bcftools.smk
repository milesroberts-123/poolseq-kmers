rule bcftools_remove_ref:
    input:
        "slim_results/{ID}.vcf",
    output:
        vcf=temp("slim_results/{ID}.vcf.gz"),
        tbi=temp("slim_results/{ID}.vcf.gz.tbi"),
        samplevcf=temp("slim_results/samples_{ID}.vcf.gz")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # compress and index
        bgzip {input}
        tabix {input}.gz
        # remove reference
        bcftools view --samples ^i0 -Oz -o {output.samplevcf} {output.vcf}
        """

rule bcftools_get_samples:
    input:
        samplevcf="slim_results/samples_{SID}.vcf.gz",
    output:
        popvcf=temp("slim_results/samples_{SID}_{PID}.vcf.gz"),
        tbi=temp("slim_results/samples_{SID}_{PID}.vcf.gz.tbi"),
        filledvcf=temp("slim_results/filled_{SID}_{PID}.vcf.gz"),
        allelefreq="slim_results/allele_freqs_{SID}_{PID}.txt",
    conda:
        "../envs/bcftools.yaml"
    params:
        sammies=get_samples,
    shell:
        """
        # get a subpopulation
        bcftools view --samples {params.sammies} -Oz -o {output.popvcf} {input}
        tabix {output.popvcf}
        # calculate allele frequencies
        bcftools +fill-tags {output.popvcf} -Oz -o {output.filledvcf}
        bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' -o {output.allelefreq} {output.filledvcf}
        """

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
