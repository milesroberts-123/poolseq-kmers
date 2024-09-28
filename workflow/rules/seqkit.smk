rule seqkit:
	input:
		vcffilled= "samples_filled_{ID}.vcf.gz",
		slimfasta = "slim_{ID}.fasta"
	output:
		samplefasta = "samples_{ID}.fasta",
		reffasta = "ref_{ID}.fasta",
		poskey = "center_kmer_pairs_{ID}.txt",
		snppos = "snp_positions_{ID}.txt"
	shell:
		"""
		seqkit grep -v -p 1 -p 2 {input.slimfasta} > {output.samplefasta}

		seqkit grep -p 1 {input.slimfasta} > {output.reffasta}

		zcat {input.vcffilled} | grep -v "^#" | cut -f 2 > {output.snppos}

		readarray posarray < {output.snppos}

		rm {output.poskey}
		
		for i in "${posarray[@]}"
		do
			echo Extracting k-mers for snp $i

			center_start=$(($i-15))
			center_end=$(($i+15))

			cat {input.slimfasta} | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> {output.poskey}
		done
		"""
