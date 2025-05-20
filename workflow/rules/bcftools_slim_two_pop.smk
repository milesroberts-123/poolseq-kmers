rule bcftools_slim_two_pop:
    input:
        compvcf = "slim_results/{ID}.vcf.gz",
        vcfidx = "slim_results/{ID}.vcf.gz.tbi",
    output:
        samplevcf_p1p2 = temp("slim_results/samples_{ID}_p1p2.vcf.gz"),
        samplevcf_p1 = temp("slim_results/samples_{ID}_p1.vcf.gz"),
        samplevcf_p2 = temp("slim_results/samples_{ID}_p2.vcf.gz"),
        allelefreq_p1 = "slim_results/allele_freqs_{ID}_p1.txt",
        allelefreq_p2 = "slim_results/allele_freqs_{ID}_p2.txt",
        filledvcf_p1 = temp("slim_results/samples_filled_{ID}_p1.vcf.gz"),
        filledvcf_p2 = temp("slim_results/samples_filled_{ID}_p2.vcf.gz"),
        filledvcf_p1p2 = temp("slim_results/samples_filled_{ID}_p1p2.vcf.gz")
    threads: 1
    resources:
        mem_mb_per_cpu=8000,
        time=239
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/bcftools_slim/{ID}.log"
    params:
        samples=get_samples,
    shell:
        """
        # remove reference
        # bcftools view --samples-file ^../config/ref.txt -Oz -o {output.samplevcf_p1p2} {input.compvcf} &>> {log}
        bcftools view --samples ^i0 -Oz -o {output.samplevcf_p1p2} {input.compvcf} &>> {log}
        
        # split by population
        bcftools view -s {params.samples} -Oz -o {output.samplevcf_p1} {output.samplevcf_p1p2} &>> {log}
        bcftools view -s ^{params.samples} -Oz -o {output.samplevcf_p2} {output.samplevcf_p1p2} &>> {log}

        # calculate allele frequencies
        bcftools +fill-tags {input.compvcf} -Oz -o {output.filledvcf_p1p2} &>> {log}
        bcftools +fill-tags {output.samplevcf_p1} -Oz -o {output.filledvcf_p1} &>> {log}
        bcftools +fill-tags {output.samplevcf_p2} -Oz -o {output.filledvcf_p2} &>> {log}

        # split allele frequencies by populations
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/MT %NS %AF %AC\n' -o {output.allelefreq_p1} {output.filledvcf_p1} &>> {log}
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/MT %NS %AF %AC\n' -o {output.allelefreq_p2} {output.filledvcf_p2} &>> {log}
        """
