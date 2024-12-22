rule varscan_two_pop:
	input:
		reffasta = "ref_{ID}_p1.fasta",
		trimbam_p2 = "trimmed_{ID}_p1.bam",
		trimbam_p1 = "trimmed_{ID}_p2.bam"
	output:
		cp1 = "calls_{ID}_p1.tsv",
		cp2 = "calls_{ID}_p2.tsv"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/varscan.yaml"
	log: 
		"logs/varscan/{ID}.log"
	shell:
		"""
		samtools mpileup -f {input.reffasta} {input.trimbam_p1} | varscan pileup2snp 1> {output.cp1} 2> {log}

		samtools mpileup -f {input.reffasta} {input.trimbam_p2} | varscan pileup2snp 1> {output.cp2} 2> {log}
		"""
