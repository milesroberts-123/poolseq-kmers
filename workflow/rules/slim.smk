def get_sigma(wildcards):
        sigma = parameters.loc[parameters["ID"] == wildcards.ID, "sigma"]
        return float(sigma.iloc[0])

def get_N(wildcards):
        N = parameters.loc[parameters["ID"] == wildcards.ID, "N"]
        return float(N.iloc[0])

def get_mu(wildcards):
        mu = parameters.loc[parameters["ID"] == wildcards.ID, "mu"]
        return float(mu.iloc[0])

def get_R(wildcards):
        R = parameters.loc[parameters["ID"] == wildcards.ID, "R"]
        return float(R.iloc[0])

def get_kappa(wildcards):
        kappa = parameters.loc[parameters["ID"] == wildcards.ID, "kappa"]
        return float(kappa.iloc[0])

rule slim:
	input:
		"../config/parameters.tsv"
	output:
		tmpVCF = temp("slim_{ID}.vcf"),
	log:
		"logs/slim/{ID}.log"
	params:
		sigma=get_sigma,
		N=get_N,
		mu=get_mu,
		R=get_R,
		kappa=get_kappa,
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/slim.yml"
	shell:
		"""
		# run simulation
		slim -d ID={wildcards.ID} -d sigma={params.sigma} -d N={params.N} -d mu={params.mu} -d R={params.R} -d kappa={params.kappa} scripts/neutral.slim &> {log}
		"""