rule poolsnp_one_pop:
	input:
		reffasta = "seqkit_results/ref_{ID}.fasta",
		trimbam = "bwa_results/{ID}.bam"
	output:
		vcf = temp("{ID}_poolsnp_output.vcf.gz"),
		cov = temp("{ID}_poolsnp_output-cov-0.9999.txt"),
		bs = temp("{ID}_poolsnp_output_BS.txt.gz"),
		mpileup = temp("{ID}.mpileup")
	params:
		wd = get_wd,
		mincount = config["mincount"]
	conda:
		"../envs/poolsnp.yaml"
	resources:
		mem_mb_per_cpu=8000,
		time=239
	benchmark:
		"benchmarks/poolsnp/{ID}.bench"
	log:
		"logs/poolsnp/{ID}.log"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam} > {output.mpileup}

		PoolSNP.sh   \
		mpileup={params.wd}{output.mpileup} \
		reference={params.wd}{input.reffasta} \
		names={wildcards.ID} \
		max-cov=0.9999 \
		min-cov={params.mincount} \
		min-count={params.mincount} \
		min-freq=0.01 \
		miss-frac=0 \
		badsites=1 \
		allsites=0 \
		output={params.wd}{wildcards.ID}_poolsnp_output &> {log}
		"""
