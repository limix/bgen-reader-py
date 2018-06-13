from bgen_reader import read_bgen

if __name__ == "__main__":

    bgen = read_bgen("complex.bgen", verbose=False)

    print(bgen["variants"].head())
    print(bgen["samples"].head())

    # Print the estimated probabilities for the first variant
    # and second individual.
    print(bgen["genotype"][0, 1].compute())

    # The NaN elements are a by-product of the heterogenous
    # ploidy and number of alleles across variants and samples.
    # For example, the 9th variant for the 4th individual
    # has ploidy
    print(bgen["X"][8, 3].compute().sel(data="ploidy").item())
    # and number of alleles equal to
    print(bgen["variants"].loc[8, "nalleles"])
    # Its probability distribution is given by the array
    print(bgen["genotype"][8, 3].compute())
    # of size
    print(len(bgen["genotype"][8, 3].compute()))

    # Since the 9th variant for the 4th individual is
    # unphased,
    print(bgen["X"][8, 3].compute().sel(data="phased").item())
    # the estimated probabilities imply the genotype

    # # Is it a phased one?
    # print(bgen["X"][0, 1].compute().sel(data="phased").item())

    # # How many haplotypes?
    # print(bgen["X"][0, 1].compute().sel(data="ploidy").item())

    # # And how many alleles?
    # print(bgen["variants"].loc[0, "nalleles"])
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # (2, 0, 0, 0, 0, 0, 0, 0)
    # # Therefore, the allele
    # print(bgen["variants"].loc[0, "allele_ids"].split(",")[1])
    # # has probability 100% of pertaining to the first haplotype.

    # # And the allele
    # print(bgen["variants"].loc[0, "allele_ids"].split(",")[0])
    # # has probability 100% of pertaining to the second haplotype.
