rule seqkit:
	input:
		vcffilled= "samples_filled_{ID}.vcf.gz",
		slimfasta = "slim_{ID}.fasta"
	output:
		samplefasta = "samples_{ID}.fasta",
		reffasta = "ref_{ID}.fasta",
		poskey = "center_kmer_pairs_{ID}.txt",
		snppos = "snp_positions_{ID}.txt"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/seqkit.yaml"
	log: 
		"logs/seqkit_{ID}.log"
	shell:
		"""
		# create file with reference genome removed
		seqkit grep -v -p 1 -p 2 {input.slimfasta} > {output.samplefasta} &> {log}

		# get a haploid reference genome
		seqkit grep -p 1 {input.slimfasta} > {output.reffasta} &> {log}
		
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
			echo Extracting k-mers for snp $i &> {log}

			center_start=$(($i-15))
			center_end=$(($i+15))

			cat {input.slimfasta} | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> {output.poskey}
		done
		"""
