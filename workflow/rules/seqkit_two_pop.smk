def get_L(wildcards):
        L = parameters.loc[parameters["ID"] == wildcards.ID, "L"]
        return int(L.iloc[0])

def get_num_genos(wildcards):
        n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
        return 2*int(n.iloc[0])

rule seqkit_two_pop:
	input:
		vcffilled= "samples_filled_{ID}_p1p2.vcf.gz",
		slimfasta = "slim_{ID}.fasta"
	output:
		tempsamplefasta = temp("samples_{ID}_p1p2.fasta"),
		p1 = "samples_{ID}_p1.fasta",
		p2 = "samples_{ID}_p2.fasta",
		reffasta = "ref_{ID}_p1.fasta",
		poskey = "center_kmer_pairs_{ID}_p1p2.txt",
		snppos = temp("snp_positions_{ID}_p1p2.txt")
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
		seqkit grep -v -n -p 1 -p 2 {input.slimfasta} > {output.tempsamplefasta}

		# get a haploid reference genome
		seqkit grep -n -p 1 {input.slimfasta} > {output.reffasta}

		# get first n individuals (population 1)
		seqkit head -n {params.n} {output.tempsamplefasta} > {output.p1}

		# get last n individuals (population 2)
		seqkit range -r -{params.n}:-1 {output.tempsamplefasta} > {output.p2}
		
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
			cat {output.tempsamplefasta} | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> {output.poskey}
		done
		"""
