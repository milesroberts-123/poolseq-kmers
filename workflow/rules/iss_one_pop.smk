def get_simtype(wildcards):
        simtype = parameters.loc[parameters["ID"] == wildcards.ID, "simtype"]
        return simtype.iloc[0]

def get_L(wildcards):
        L = parameters.loc[parameters["ID"] == wildcards.ID, "L"]
        return int(L.iloc[0])

def get_cov(wildcards):
        cov = parameters.loc[parameters["ID"] == wildcards.ID, "cov"]
        return int(cov.iloc[0])

def get_sequencer(wildcards):
        sequencer = parameters.loc[parameters["ID"] == wildcards.ID, "sequencer"]
        return sequencer.iloc[0]

rule iss_one_pop:
	input:
		"samples_{ID}.fasta"
		#p1="samples_{ID}_p1.fasta",
		#p2="samples_{ID}_p2.fasta"
	output:
		#temp("reads_{ID}_p1_R1.fastq"),
		#temp("reads_{ID}_p1_R2.fastq"),
		#temp("reads_{ID}_p2_R1.fastq"),
		#temp("reads_{ID}_p2_R2.fastq")
		temp("reads_{ID}_R1.fastq"),
		temp("reads_{ID}_R2.fastq")
	threads: 4
	resources:
		mem_mb_per_cpu=2000,
		time=239
	conda:
		"../envs/iss.yaml"
	log:
		"logs/iss/{ID}.log"
	params:
		L = get_L,
		cov = get_cov,
		sequencer = get_sequencer
	shell:
		"""
		# read length = 300 bp
		if [ "{params.sequencer}" == "miseq" ] || [ "{params.sequencer}" == "nextseq" ]; then
			# calculate number of reads for desired coverage level
			nreads=$(({params.L}*{params.cov}/300))
		fi

		# read length 150 bp
		if [ "{params.sequencer}" == "novaseq" ]; then
			# calculate number of reads for desired coverage level
			nreads=$(({params.L}*{params.cov}/150))
		fi
		
		# read length 125 bp
		if [ "{params.sequencer}" == "hiseq" ]; then
			# calculate number of reads for desired coverage level
			nreads=$(({params.L}*{params.cov}/125))
		fi

		# simulate reads
		iss generate -g {input} --cpus {threads} --model {params.sequencer} -n $nreads --abundance uniform --output reads_{wildcards.ID} &> {log}
		"""
