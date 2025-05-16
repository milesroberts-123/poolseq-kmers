def get_simtype(wildcards):
        simtype = parameters.loc[parameters["ID"] == wildcards.ID, "simtype"]
        return simtype.iloc[0]

def get_sigma(wildcards):
        sigma = parameters.loc[parameters["ID"] == wildcards.ID, "sigma"]
        return float(sigma.iloc[0])

def get_N(wildcards):
        N = parameters.loc[parameters["ID"] == wildcards.ID, "N"]
        return int(N.iloc[0])

def get_N1(wildcards):
        N1 = parameters.loc[parameters["ID"] == wildcards.ID, "N1"]
        return int(N1.iloc[0])

def get_N2(wildcards):
        N2 = parameters.loc[parameters["ID"] == wildcards.ID, "N2"]
        return int(N2.iloc[0])

def get_mg1(wildcards):
        mg1 = parameters.loc[parameters["ID"] == wildcards.ID, "mg1"]
        return float(mg1.iloc[0])

def get_mg2(wildcards):
        mg2 = parameters.loc[parameters["ID"] == wildcards.ID, "mg2"]
        return float(mg2.iloc[0])

def get_h(wildcards):
        h = parameters.loc[parameters["ID"] == wildcards.ID, "h"]
        return float(h.iloc[0])

def get_s(wildcards):
        s = parameters.loc[parameters["ID"] == wildcards.ID, "s"]
        return float(s.iloc[0])

def get_n(wildcards):
        n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
        return int(n.iloc[0])

def get_mu(wildcards):
        mu = parameters.loc[parameters["ID"] == wildcards.ID, "mu"]
        return float(mu.iloc[0])

def get_R(wildcards):
        R = parameters.loc[parameters["ID"] == wildcards.ID, "R"]
        return float(R.iloc[0])

#def get_L(wildcards):
#        L = parameters.loc[parameters["ID"] == wildcards.ID, "L"]
#        return int(L.iloc[0])

def get_tau(wildcards):
        tau = parameters.loc[parameters["ID"] == wildcards.ID, "tau"]
        return int(tau.iloc[0])

def get_qtl_mean(wildcards):
        qtl_mean = parameters.loc[parameters["ID"] == wildcards.ID, "qtl_mean"]
        return float(qtl_mean.iloc[0])

def get_qtl_sigma(wildcards):
        qtl_sigma = parameters.loc[parameters["ID"] == wildcards.ID, "qtl_sigma"]
        return float(qtl_sigma.iloc[0])

def get_qtl_prop(wildcards):
        qtl_prop = parameters.loc[parameters["ID"] == wildcards.ID, "qtl_prop"]
        return float(qtl_prop.iloc[0])

def get_optimum_mean(wildcards):
        optimum_mean = parameters.loc[parameters["ID"] == wildcards.ID, "optimum_mean"]
        return float(optimum_mean.iloc[0])

def get_optimum_sigma(wildcards):
        optimum_sigma = parameters.loc[parameters["ID"] == wildcards.ID, "optimum_sigma"]
        return float(optimum_sigma.iloc[0])

def get_phenotype_cutoff(wildcards):
        phenotype_cutoff = parameters.loc[parameters["ID"] == wildcards.ID, "phenotype_cutoff"]
        return float(phenotype_cutoff.iloc[0])

rule slim:
	input:
		"ancestral_genome_results/{ID}.fasta"
	output:
		temp("slim_results/{ID}.vcf"),
		temp("slim_results/{ID}.fasta")
	log:
		"logs/slim/{ID}.log"
	params:
		simtype=get_simtype,
		sigma=get_sigma,
		N=get_N,
		N1=get_N1,
		N2=get_N2,
		mg1=get_mg1,
		mg2=get_mg2,
		n=get_n,
		h=get_h,
		s=get_s,
		mu=get_mu,
		R=get_R,
		#L=get_L,
		tau=get_tau,
		qtl_mean=get_qtl_mean,
		qtl_sigma=get_qtl_sigma,
		qtl_prop=get_qtl_prop,
		optimum_mean=get_optimum_mean,
		optimum_sigma=get_optimum_sigma,
		phenotype_cutoff=get_phenotype_cutoff
	threads: 1
	resources:
		mem_mb_per_cpu=8000,
		time=239
	conda:
		"../envs/slim.yaml"
	shell:
		"""
		if [ "{params.simtype}" == "onepop" ]; then
			slim -d ID={wildcards.ID} -d sigma={params.sigma} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/neutral.slim &> {log}
		fi

		if [ "{params.simtype}" == "twopop" ]; then
			slim -d ID={wildcards.ID} -d sigma={params.sigma} -d N1={params.N1} -d N2={params.N2} -d mg1={params.mg1} -d mg2={params.mg2} -d mu={params.mu} -d R={params.R} -d n={params.n} -d tau={params.tau} scripts/two_pop.slim &> {log}
		fi

		if [ "{params.simtype}" == "bsa" ]; then
			slim -d ID={wildcards.ID} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} -d qtl_mean={params.qtl_mean} -d qtl_sigma={params.qtl_sigma} -d qtl_prop={params.qtl_prop} -d optimum_mean={params.optimum_mean} -d optimum_sigma={params.optimum_sigma} -d phenotype_cutoff={params.phenotype_cutoff} scripts/bsa.slim &> {log}
		fi

		if [ "{params.simtype}" == "sweep" ]; then
			slim -d ID={wildcards.ID} -d h={params.h} -d s={params.s} -d sigma={params.sigma} -d N={params.N} -d mu={params.mu} -d R={params.R} -d n={params.n} scripts/sweep.slim &> {log}
		fi
		"""
