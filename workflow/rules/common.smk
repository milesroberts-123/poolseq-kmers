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

def get_shuffle(wildcards):
    shuffle = parameters.loc[parameters["ID"] == wildcards.ID, "shuffle"]
    return shuffle.iloc[0]

def get_samples(wildcards):
    # get sample size
    n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
    n = int(n.iloc[0])

    # create list of sample names from slim convention
    samples = list(range(1, n + 1))
    samples = ["i" + str(x) for x in samples] 

    # create comma-sep list for bcftools
    samples = ','.join(samples)

    return samples

# n is number of individuals
# multiply by 2 to convert to number of genomes
def get_pool(wildcards):
    n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
    return 2*int(n.iloc[0])

def get_cov(wildcards):
    cov = parameters.loc[parameters["ID"] == wildcards.ID, "cov"]
    return int(cov.iloc[0])

def get_sequencer(wildcards):
    sequencer = parameters.loc[parameters["ID"] == wildcards.ID, "sequencer"]
    return sequencer.iloc[0]

def get_wd(wildcards):
    return os.getcwd() + "/"

def get_num_genos(wildcards):
    n = parameters.loc[parameters["ID"] == wildcards.ID, "n"]
    return int(n.iloc[0])
