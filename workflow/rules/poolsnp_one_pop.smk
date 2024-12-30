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

rule poolsnp_one_pop:
	input:
		reffasta = "ref_{ID}.fasta",
		trimbam = "trimmed_{ID}.bam"
	output:
		vcf = "{ID}_poolsnp_output.vcf.gz",
		cov = "{ID}_poolsnp_output-cov-0.98.txt",
		bs = "{ID}_poolsnp_output_BS.txt.gz",
		mpileup = temp("{ID}.mpileup")
	params:
		#names = get_names,
		wd = get_wd,
		prefix = get_prefix
	conda:
		"../envs/poolsnp.yaml"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam} > {output.mpileup}

		PoolSNP.sh   \
		mpileup={params.wd}{output.mpileup} \
		reference={params.wd}{input.reffasta} \
		names={wildcards.ID} \
		max-cov=0.98 \
		min-cov=10 \
		min-count=10 \
		min-freq=0.01 \
		miss-frac=0.2 \
		badsites=1 \
		allsites=0 \
		output={params.prefix}
		"""