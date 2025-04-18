def get_pA(wildcards):
        pA = parameters.loc[parameters["ID"] == wildcards.ID, "pA"]
        return float(pA.iloc[0])

def get_pC(wildcards):
        pC = parameters.loc[parameters["ID"] == wildcards.ID, "pC"]
        return float(pC.iloc[0])

def get_pG(wildcards):
        pG = parameters.loc[parameters["ID"] == wildcards.ID, "pG"]
        return float(pG.iloc[0])

def get_pT(wildcards):
        pT = parameters.loc[parameters["ID"] == wildcards.ID, "pT"]
        return float(pT.iloc[0])

def get_shape(wildcards):
        shape = parameters.loc[parameters["ID"] == wildcards.ID, "shape"]
        return float(shape.iloc[0])

def get_L(wildcards):
        L = parameters.loc[parameters["ID"] == wildcards.ID, "L"]
        return int(L.iloc[0])

rule ancestral_genome:
	input:
		"../config/parameters.tsv"
	output:
		"ancestral_seq_{ID}.fasta",
		#"power_law_{ID}.jpg"
	log:
		"logs/ancestral_genome/{ID}.log"
	params:
		pA=get_pA,
		pC=get_pC,
		pG=get_pG,
		pT=get_pT,
		shape=get_shape,
		L=get_L,
	conda:
		"../envs/R.yaml"
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	shell:
		"Rscript scripts/random_genome.R {params.pA} {params.pC} {params.pG} {params.pT} {params.shape} 31 {params.L} {wildcards.ID} &> {log}"