rule blast:
	input:
		ref = "seqkit_results/ref_{ID}.fasta",
		unitigs = "unitig_caller_results/unitigs_renamed_{ID}.fasta"
	output:
		ndb = temp("seqkit_results/ref_{ID}.ndb"),
		nhr = temp("seqkit_results/ref_{ID}.nhr"),
		nin = temp("seqkit_results/ref_{ID}.nin"),
		njs = temp("seqkit_results/ref_{ID}.njs"),
		not_dbfile = temp("seqkit_results/ref_{ID}.not"),
		nsq = temp("seqkit_results/ref_{ID}.nsq"),
		ntf = temp("seqkit_results/ref_{ID}.ntf"),
		nto = temp("seqkit_results/ref_{ID}.nto"),
		alignments = "blast_results/{ID}.txt"
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