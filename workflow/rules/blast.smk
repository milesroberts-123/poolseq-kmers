rule blast:
	input:
		ref = "ref_{ID}.fasta",
		unitigs = "unitigs_renamed_{ID}.fasta"
	output:
		ndb = temp("ref_{ID}.ndb"),
		nhr = temp("ref_{ID}.nhr"),
		nin = temp("ref_{ID}.nin"),
		njs = temp("ref_{ID}.njs"),
		not_dbfile = temp("ref_{ID}.not"),
		nsq = temp("ref_{ID}.nsq"),
		ntf = temp("ref_{ID}.ntf"),
		nto = temp("ref_{ID}.nto"),
		alignments = "unitig_alignments_{ID}.txt"
	conda:
		"../envs/blast.yaml"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	shell:
		"""
		makeblastdb -in {input.ref} -title $(basename {input.ref} .fasta) -dbtype nucl -out $(basename {input.ref} .fasta)

		blastn -query {input.unitigs} -db $(basename {input.ref} .fasta) -out {output.alignments} -max_target_seqs 1 -evalue 1e-10 -outfmt 6
		"""