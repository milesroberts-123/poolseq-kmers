rule ancestral_genome:
    input:
        "../config/parameters.tsv",
    output:
        temp("ancestral_genome_results/{ID}.fasta"),
    params:
        pA=lookup(query="ID == '{ID}'", within=parameters, cols="pA"),
        pC=lookup(query="ID == '{ID}'", within=parameters, cols="pC"),
        pG=lookup(query="ID == '{ID}'", within=parameters, cols="pG"),
        pT=lookup(query="ID == '{ID}'", within=parameters, cols="pT"),
        shape=lookup(query="ID == '{ID}'", within=parameters, cols="shape"),
        L=lookup(query="ID == '{ID}'", within=parameters, cols="L"),
        k=lookup(query="ID == '{ID}'", within=parameters, cols="k"),
        shuffleKmers=lookup(query="ID == '{ID}'", within=parameters, cols="shuffle"),
    conda:
        "../envs/R.yaml"
    shell:
        """
        Rscript scripts/random_genome.R {params.pA} {params.pC} {params.pG} {params.pT} {params.shape} {params.k} {params.L} {wildcards.ID} {params.shuffleKmers}
        """

rule exact_kmc_histo:
    input:
        "ancestral_genome_results/{ID}.fasta",
    output:
        histo="exact_kmc_histo_results/{ID}.histo",
        pre=temp("counts_{ID}.kmc_pre"),
        suf=temp("counts_{ID}.kmc_suf"),
    conda:
        "../envs/kmc.yaml"
    params:
        k=lookup(query="ID == '{ID}'", within=parameters, cols="k"),
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}" ]; then
            rm -r tmp_kmc_{wildcards.ID}
        fi

        mkdir tmp_kmc_{wildcards.ID}

        # count k-mers
        kmc -t{threads} -m9 -ci1 -cs100000 -fm -k{params.k} {input} counts_{wildcards.ID} tmp_kmc_{wildcards.ID}

        # convert to histogram
        kmc_tools transform counts_{wildcards.ID} histogram {output.histo}

        # rm tmp dir
        rm -r tmp_kmc_{wildcards.ID}
        """

rule slim:
    input:
        "ancestral_genome_results/{ID}.fasta",
    output:
        temp("slim_results/{ID}.vcf"),
        temp("slim_results/{ID}.fasta"),
    params:
        simtype=lookup(query="ID == '{ID}'", within=parameters, cols="simtype"),
        slimseed=lookup(query="ID == '{ID}'", within=parameters, cols="slimseed"),
        N=lookup(query="ID == '{ID}'", within=parameters, cols="N"),
        n=lookup(query="ID == '{ID}'", within=parameters, cols="n"),
        h=lookup(query="ID == '{ID}'", within=parameters, cols="h"),
        s=lookup(query="ID == '{ID}'", within=parameters, cols="s"),
        mu=lookup(query="ID == '{ID}'", within=parameters, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=parameters, cols="R"),
        N1=lookup(query="ID == '{ID}'", within=parameters, cols="N1"),
        N2=lookup(query="ID == '{ID}'", within=parameters, cols="N2"),
        mg1=lookup(query="ID == '{ID}'", within=parameters, cols="mg1"),
        mg2=lookup(query="ID == '{ID}'", within=parameters, cols="mg2"),
        tau=lookup(query="ID == '{ID}'", within=parameters, cols="tau"),
    conda:
        "../envs/slim.yaml"
    shell:
        """
        if [ "{params.simtype}" == "onepop" ]; then
            slim -d ID={wildcards.ID} -d SLIMSEED={params.slimseed} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/neutral.slim
        fi

        if [ "{params.simtype}" == "sweep" ]; then
            slim -d ID={wildcards.ID} -d SLIMSEED={params.slimseed} -d h={params.h} -d s={params.s} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/sweep.slim
        fi

        if [ "{params.simtype}" == "twopop" ]; then
            slim -d ID={wildcards.ID} -d SLIMSEED={params.slimseed} -d N1={params.N1} -d N2={params.N2} -d mg1={params.mg1} -d mg2={params.mg2} -d mu={params.mu} -d R={params.R} -d n={params.n} -d tau={params.tau} scripts/two_pop.slim
        fi
        """

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

def sample_start(wildcards):
    # get sample size
    n = parameters.loc[parameters["ID"] == wildcards.SID, "n"]
    n = int(n.iloc[0])

    # get population
    p = int(wildcards.PID)

    return 1 + n * p


def sample_end(wildcards):
    # get sample size
    n = parameters.loc[parameters["ID"] == wildcards.SID, "n"]
    n = int(n.iloc[0])

    # get population
    p = int(wildcards.PID)

    return n * (p + 1)


rule seqkit_get_samples:
    input:
        slimfasta="slim_results/{SID}.fasta",
    output:
        tempsamplefasta=temp("seqkit_results/temp_{SID}_{PID}.fasta"),
        pop1="seqkit_results/samples_{SID}_{PID}.fasta",
    params:
        start=sample_start,
        end=sample_end,
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # create file with reference genome removed
        seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.tempsamplefasta}

        # get group of n individuals
        seqkit range -r {params.start}:{params.end} {output.tempsamplefasta} > {output.pop1}        
        """

def calc_nreads(wildcards):
    L = parameters.loc[parameters["ID"] == wildcards.SID, "L"]
    cov = parameters.loc[parameters["ID"] == wildcards.SID, "cov"]
    sequencer = parameters.loc[parameters["ID"] == wildcards.SID, "sequencer"]

    L = int(L.iloc[0])
    cov = int(cov.iloc[0])
    sequencer = sequencer.iloc[0]

    if sequencer == "miseq" or sequencer == "nextseq":
        nreads=(L*cov)/300
    if sequencer == "novaseq":
        nreads=(L*cov)/150
    if sequencer == "hiseq":
        nreads=(L*cov)/125

    return int(nreads)

rule iss:
    input:
        "seqkit_results/samples_{SID}_{PID}.fasta",
    output:
        temp("iss_results/reads_{SID}_{PID}_R1.fastq"),
        temp("iss_results/reads_{SID}_{PID}_R2.fastq"),
    conda:
        "../envs/iss.yaml"
    params:
        sequencer=lookup(query="ID == '{SID}'", within=parameters, cols="sequencer"),
        issseed=lookup(query="ID == '{SID}'", within=parameters, cols="issseed"),
        nreads=calc_nreads
    shell:
        """
        echo Number of reads to simulate: 
        echo {params.nreads}

        # simulate reads
        iss generate -g {input} --seed {params.issseed} --cpus {threads} --model {params.sequencer} -n {params.nreads} --abundance uniform --output iss_results/reads_{wildcards.SID}_{wildcards.PID}
        """

rule fastp:
    input:
        read1="iss_results/reads_{SID}_{PID}_R1.fastq",
        read2="iss_results/reads_{SID}_{PID}_R2.fastq",
    output:
        pread1=temp("fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq"),
        pread2=temp("fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{SID}_{PID}.fastq"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{SID}_{PID}.fastq"),
        jsonR1R2="fastp_results/{SID}_{PID}_R1R2.json",
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=lookup(query="ID == '{SID}'", within=parameters, cols="k"),
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        # trim low quality bases
        fastp --thread {threads} -u {params.unqualLimit} -q {params.qualThresh} --correction -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {input.read1} -I {input.read2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2}
        """

rule ancestral_samtools_faidx:
    input:
        "ancestral_genome_results/{SID}.fasta",
    output:
        temp("ancestral_genome_results/{SID}.fasta.fai"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # index reference genome
        samtools faidx {input}
        """


rule freqk_index:
    input:
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        fasta="ancestral_genome_results/{SID}.fasta",
        fai="ancestral_genome_results/{SID}.fasta.fai",
    output:
        index=temp("freqk_indices/{SID}_{PID}.txt"),
    benchmark:
        "benchmarks/freqk_index/{SID}_{PID}.bench"
    params:
        k=lookup(query="ID == '{SID}'", within=parameters, cols="k"),
    shell:
        """
        # index panel of variants
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """


rule freqk_var_dedup:
    input:
        "freqk_indices/{SID}_{PID}.txt",
    output:
        temp("freqk_var_dedup/{SID}_{PID}.txt"),
    benchmark:
        "benchmarks/freqk_var_dedup/{SID}_{PID}.bench"
    shell:
        """
        ./scripts/freqk var-dedup --index {input} --output {output}
        """


rule freqk_ref_dedup:
    input:
        index="freqk_var_dedup/{SID}_{PID}.txt",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        fasta="ancestral_genome_results/{SID}.fasta",
        fai="ancestral_genome_results/{SID}.fasta.fai",
    output:
        "freqk_ref_dedup/{SID}_{PID}.txt",
    benchmark:
        "benchmarks/freqk_ref_dedup/{SID}_{PID}.bench"
    shell:
        """
        ./scripts/freqk ref-dedup --index {input.index} --fasta {input.fasta} --vcf {input.vcf} --output {output}
        """


rule combine_fastqs:
    input:
        pread1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
        uread1="fastp_results/trimmed_unpaired_R1_{SID}_{PID}.fastq",
        uread2="fastp_results/trimmed_unpaired_R2_{SID}_{PID}.fastq",
    output:
        temp("all_{SID}_{PID}.fastq"),
    shell:
        "cat {input.pread1} {input.pread2} {input.uread1} {input.uread2} > {output}"


rule freqk_count:
    input:
        reads="all_{SID}_{PID}.fastq",
        index="freqk_ref_dedup/{SID}_{PID}.txt",
    output:
        counts="freqk_results/{SID}_{PID}_counts.txt",
        freqs="freqk_results/{SID}_{PID}_freqs.txt",
    benchmark:
        "benchmarks/freqk_count/{SID}_{PID}.bench"
    shell:
        """
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {input.reads} --freq-output {output.freqs} --count-output {output.counts}
        """

rule vg_autoindex:
    input:
        fasta="ancestral_genome_results/{SID}.fasta",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
    output:
        dist=temp("{SID}_{PID}.dist"),
        gbz=temp("{SID}_{PID}.giraffe.gbz"),
        min=temp("{SID}_{PID}.shortread.withzip.min"),
        zip=temp("{SID}_{PID}.shortread.zipcodes"),
    benchmark:
        "benchmarks/vg_autoindex/{SID}_{PID}.bench"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe --threads {threads} -r {input.fasta} -v {input.vcf} -p {wildcards.SID}_{wildcards.PID}"


rule vg_giraffe_paired:
    input:
        dist="{SID}_{PID}.dist",
        gbz="{SID}_{PID}.giraffe.gbz",
        min="{SID}_{PID}.shortread.withzip.min",
        zip="{SID}_{PID}.shortread.zipcodes",
        pread1="fastp_results/trimmed_paired_R1_{SID}_{PID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/paired_{SID}_{PID}.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{SID}_{PID}.bench"
    shell:
        "vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -t {threads} -f {input.pread1} -f {input.pread2} > {output}"


rule vg_giraffe_unpaired:
    input:
        dist="{SID}_{PID}.dist",
        gbz="{SID}_{PID}.giraffe.gbz",
        min="{SID}_{PID}.shortread.withzip.min",
        zip="{SID}_{PID}.shortread.zipcodes",
        uread="fastp_results/trimmed_unpaired_{read}_{SID}_{PID}.fastq",
    output:
        temp("vg_giraffe_results/unpaired_{read}_{SID}_{PID}.gam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/unpaired_{read}_{SID}_{PID}.bench"
    shell:
        "vg giraffe -Z {input.gbz} -d {input.dist} -m {input.min} -z {input.zip} -t {threads} -f {input.uread} > {output}"


rule vg_surject_unpaired:
    input:
        gbz="{SID}_{PID}.giraffe.gbz",
        gam="vg_giraffe_results/unpaired_{read}_{SID}_{PID}.gam",
    output:
        temp("vg_surject_results/unpaired_{read}_{SID}_{PID}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_surject/unpaired_{read}_{SID}_{PID}.bench"
    shell:
        "vg surject -x {input.gbz} -t {threads} -b {input.gam} > {output}"


rule vg_surject_paired:
    input:
        gbz="{SID}_{PID}.giraffe.gbz",
        gam="vg_giraffe_results/paired_{SID}_{PID}.gam",
    output:
        temp("vg_surject_results/paired_{SID}_{PID}.bam"),
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_surject/paired_{SID}_{PID}.bench"
    shell:
        "vg surject -x {input.gbz} -t {threads} -b {input.gam} > {output}"


rule samtools_sort_unpaired:
    input:
        "vg_surject_results/unpaired_{read}_{SID}_{PID}.bam",
    output:
        temp("samtools_sort_results/unpaired_{read}_{SID}_{PID}.bam"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_sort/unpaired_{read}_{SID}_{PID}.bench"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_sort_paired:
    input:
        "vg_surject_results/paired_{SID}_{PID}.bam",
    output:
        temp("samtools_sort_results/paired_{SID}_{PID}.bam"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_sort/paired_{SID}_{PID}.bench"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_merge:
    input:
        "samtools_sort_results/paired_{SID}_{PID}.bam",
        expand(
            "samtools_sort_results/unpaired_{read}_{{SID}}_{{PID}}.bam",
            read=["R1", "R2"],
        ),
    output:
        temp("samtools_merge_results/{SID}_{PID}.bam"),
    benchmark:
        "benchmarks/samtools_merge/{SID}_{PID}.bench"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "samtools merge -o {output} {input}"

rule samtools_markdup:
    input:
        "samtools_merge_results/{SID}_{PID}.bam"
    output:
        temp("samtools_markdup_results/{SID}_{PID}.bam")
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/samtools_markdup/{SID}_{PID}.bench"
    shell:
        "samtools collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} - {output}"

rule freebayes_vg:
    input:
        reffasta="ancestral_genome_results/{SID}.fasta",
        index="ancestral_genome_results/{SID}.fasta.fai",
        trimbam="samtools_markdup_results/{SID}_{PID}.bam",
        vcf="slim_results/samples_{SID}_{PID}.vcf.gz",
        tbi="slim_results/samples_{SID}_{PID}.vcf.gz.tbi",
    output:
        "freebayes_vg_results/{SID}_{PID}.vcf",
    conda:
        "../envs/freebayes.yaml"
    benchmark:
        "benchmarks/freebayes_vg/{SID}_{PID}.bench"
    params:
        n=get_pool,
    shell:
        "freebayes -f {input.reffasta} -p {params.n} --use-best-n-alleles 2 --variant-input {input.vcf} --only-use-input-alleles --pooled-discrete {input.trimbam} 1> {output}"

rule varscan_vg:
    input:
        reffasta="ancestral_genome_results/{SID}.fasta",
        index="ancestral_genome_results/{SID}.fasta.fai",
        trimbam="samtools_markdup_results/{SID}_{PID}.bam",
    output:
        "varscan_results/{SID}_{PID}.tsv",
    conda:
        "../envs/varscan.yaml"
    benchmark:
        "benchmarks/varscan/{SID}_{PID}.bench"
    params:
        mincount=config["mincount"],
        minfreq=config["minfreq"]
    shell:
        "samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp --min-coverage {params.mincount} --min-var-freq {params.minfreq} 1> {output}"

rule poolsnp_vg:
    input:
        reffasta="ancestral_genome_results/{SID}.fasta",
        index="ancestral_genome_results/{SID}.fasta.fai",
        trimbam="samtools_markdup_results/{SID}_{PID}.bam",
    output:
        vcf=temp("{SID}_{PID}_poolsnp_output.vcf.gz"),
        cov=temp("{SID}_{PID}_poolsnp_output-cov-0.9999.txt"),
        bs=temp("{SID}_{PID}_poolsnp_output_BS.txt.gz"),
        mpileup=temp("{SID}_{PID}.mpileup"),
    params:
        wd=get_wd,
        mincount=config["mincount"],
        minfreq=config["minfreq"]
    conda:
        "../envs/poolsnp.yaml"
    benchmark:
        "benchmarks/poolsnp/{SID}_{PID}.bench"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} > {output.mpileup}

        PoolSNP.sh   \
        mpileup={params.wd}{output.mpileup} \
        reference={params.wd}{input.reffasta} \
        names={wildcards.SID} \
        max-cov=0.9999 \
        min-cov={params.mincount} \
        min-count=1 \
        min-freq={params.minfreq} \
        miss-frac=0 \
        badsites=1 \
        allsites=0 \
        output={params.wd}{wildcards.SID}_{wildcards.PID}_poolsnp_output
        """
