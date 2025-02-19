def get_L(wildcards):
        L = parameters.loc[parameters["ID"] == wildcards.ID, "L"]
        return int(L.iloc[0])

def get_num_genos(wildcards):
        n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
        return int(n.iloc[0])

rule seqkit_one_pop:
	input:
		vcffilled= "samples_filled_{ID}.vcf.gz",
		slimfasta = "slim_{ID}.fasta"
	output:
		samplefasta = "samples_{ID}.fasta",
		reffasta = "ref_{ID}.fasta",
		poskey = "center_kmer_pairs_{ID}.txt",
		snppos = temp("snp_positions_{ID}.txt")
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	params:
		L=get_L,
		n=get_num_genos
	conda:
		"../envs/seqkit.yaml"
	log: 
		"logs/seqkit/{ID}.log"
	shell:
		"""
		# create file with reference genome removed
		seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.samplefasta}

		# get a haploid reference genome
		seqkit grep -n -p 1 {input.slimfasta} > {output.reffasta}
		
		# get list of snp positions
		zcat {input.vcffilled} | grep -v "^#" | cut -f 2 > {output.snppos}

		# load list
		readarray posarray < {output.snppos}

		# if output already exists (i.e. previous failed run) then delete it
		if [ -e "{output.poskey}" ]; then
			rm {output.poskey}
		fi
		
		# loop over snp positions
		for i in "${{posarray[@]}}"
		do
			echo Extracting k-mers for snp $i &>> {log}

			# calculate bounds of k-mers
			center_start=$(($i-15))
			center_end=$(($i+15))

			# check if bounds of k-mer extend beyond edges of chromosome
			if [ "$center_start" -lt 1 ]; then
				continue
			fi

			if [ "$center_end" -gt {params.L} ]; then
				continue
			fi
			
			# extract k-mers
			cat {output.samplefasta} | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> {output.poskey}
		done
		"""
