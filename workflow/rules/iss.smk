def get_L(wildcards):
        L = parameters.loc[parameters["ID"] == wildcards.ID, "L"]
        return float(L.iloc[0])

def get_cov(wildcards):
        cov = parameters.loc[parameters["ID"] == wildcards.ID, "cov"]
        return float(cov.iloc[0])

rule iss:
	input:
		"samples_{ID}.fasta"
	output:
		"reads_{ID}_R1.fastq",
		"reads_{ID}_R2.fastq"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/iss.yaml"
	log:
		"logs/iss_{ID}.log"
	params:
		L = get_L,
		cov = get_cov
	shell:
		"""
		# calculate number of reads for desired coverage level
		nreads=$(({params.L}*{params.cov}/300))

		# simulate reads
		iss generate -g {input} --cpus {threads} --model miseq -n $nreads --abundance uniform --output reads_{wildcards.ID} &> {log}
		"""
