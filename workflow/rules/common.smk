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

def get_wd(wildcards):
    return os.getcwd() + "/"
