from bgen_reader import read_bgen

if __name__ == "__main__":

    bgen = read_bgen("haplotypes.bgen", verbose=False)

    print(bgen["variants"].head())
    print(bgen["samples"].head())

    # Print the estimated probabilities for the first variant
    # and second individual.
    print(bgen["genotype"][0, 1].compute())

    # Is it a phased one?
    print(bgen["X"][0, 1].compute().sel(data="phased").item())

    # How many haplotypes?
    print(bgen["X"][0, 1].compute().sel(data="ploidy").item())

    # And how many alleles?
    print(bgen["variants"].loc[0, "nalleles"])

    # Therefore, the allele
    print(bgen["variants"].loc[0, "allele_ids"].split(",")[1])
    # has probability 100% of pertaining to the first haplotype.

    # And the allele
    print(bgen["variants"].loc[0, "allele_ids"].split(",")[0])
    # has probability 100% of pertaining to the second haplotype.
