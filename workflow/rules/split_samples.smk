def get_num_genos(wildcards):
        n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
        return 2*int(n.iloc[0])

rule split_samples:
	input:
		"samples_{ID}.fasta"
	output:
		p1=temp("samples_{ID}_p1.fasta"),
		p2=temp("samples_{ID}_p2.fasta")
	log:
		"logs/split_samples/{ID}.log"
	params:
		simtype=get_simtype,
		n=get_num_genos
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/seqkit.yaml"
	wildcard_constraints:
		ID="^_"
	shell:
		"""
		seqkit head -n {params.n} {input} > {output.p1}
		seqkit range -r -{params.n}:-1 {input} > {output.p2}
		"""
	