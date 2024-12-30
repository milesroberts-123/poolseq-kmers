def get_samples(wildcards):
	# get sample size
        n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
	n = int(n.iloc[0])

	# create list of sample names from slim convention
	samples = list(range(1, n + 1))
	samples = ["i" + str(x) for x in samples] 

	# create comma-sep list for bcftools
	samples = ','.join(samples)

        return samples

def get_wd(wildcards):
	return os.getcwd() + "/"

def get_prefix(wildcards):
	return os.getcwd() + "/" + str(wildcards.ID) + "_poolsnp_output" 

rule poolsnp_two_pop:
	input:
		reffasta = "ref_{ID}_p1.fasta",
		trimbam_p1 = "trimmed_{ID}_p1.bam",
		trimbam_p2 = "trimmed_{ID}_p2.bam"
	output:
		vcf_p1 = "{ID}_p1_poolsnp_output.vcf.gz",
		cov_p1 = "{ID}_p1_poolsnp_output-cov-0.98.txt",
		bs_p1 = "{ID}_p1_poolsnp_output_BS.txt.gz",
		mpileup_p1 = temp("{ID}_p1.mpileup"),
		vcf_p2 = "{ID}_p2_poolsnp_output.vcf.gz",
		cov_p2 = "{ID}_p2_poolsnp_output-cov-0.98.txt",
		bs_p2 = "{ID}_p2_poolsnp_output_BS.txt.gz",
		mpileup_p2 = temp("{ID}_p2.mpileup")
	params:
		#names = get_names,
		wd = get_wd,
		prefix = get_prefix
	conda:
		"../envs/poolsnp.yaml"
	resources:
		mem_mb_per_cpu=8000,
		time=239
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam_p1} > {output.mpileup_p1}

		samtools mpileup -f {input.reffasta} {input.trimbam_p2} > {output.mpileup_p2}

		PoolSNP.sh mpileup={params.wd}{output.mpileup_p1} reference={params.wd}{input.reffasta} names={wildcards.ID}_p1 max-cov=0.98 min-cov=10 min-count=10 min-freq=0.01 miss-frac=0.2 badsites=1 allsites=0 output={params.prefix}_p1
	
        	PoolSNP.sh mpileup={params.wd}{output.mpileup_p2} reference={params.wd}{input.reffasta} names={wildcards.ID}_p2 max-cov=0.98 min-cov=10 min-count=10 min-freq=0.01 miss-frac=0.2 badsites=1 allsites=0 output={params.prefix}_p2
		"""
