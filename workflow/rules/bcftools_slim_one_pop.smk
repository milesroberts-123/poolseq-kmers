rule bcftools_slim_one_pop:
    input:
        compvcf = "slim_results/{ID}.vcf.gz",
        vcfidx = "slim_results/{ID}.vcf.gz.tbi",
    output:
        samplevcf = temp("slim_results/samples_{ID}.vcf.gz"),
        allelefreq = "slim_results/allele_freqs_{ID}.txt",
        filledvcf = temp("slim_results/samples_filled_{ID}.vcf.gz")
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/bcftools_slim/{ID}.log"
    shell:
        """
        # remove reference
        # bcftools view --samples-file ^../config/ref.txt -Oz -o {output.samplevcf} {input.compvcf} &>> {log}
        bcftools view --samples ^i0 -Oz -o {output.samplevcf} {input.compvcf} &>> {log}

        # calculate allele frequencies
        bcftools +fill-tags {output.samplevcf} -Oz -o {output.filledvcf} &>> {log}
        bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' -o {output.allelefreq} {output.filledvcf} &>> {log}
        """
