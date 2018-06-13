from bgen_reader import read_bgen, convert_to_dosage

if __name__ == "__main__":

    bgen = read_bgen("complex.bgen", verbose=False)

    print(bgen["variants"])
    print(bgen["samples"])

    # Print the estimated probabilities for the first variant
    # and second individual.
    print(bgen["genotype"][0, 1].compute())

    # The NaN elements are a by-product of the heterogenous
    # ploidy and number of alleles across variants and samples.
    # For example, the 9th variant for the 4th individual
    # has ploidy
    ploidy = bgen["X"][8, 3].compute().sel(data="ploidy").item()
    print(ploidy)
    # and number of alleles equal to
    nalleles = bgen["variants"].loc[8, "nalleles"]
    print(nalleles)
    # Its probability distribution is given by the array
    p = bgen["genotype"][8, 3].compute()
    print(p)
    # of size
    print(len(p))

    # Since the 9th variant for the 4th individual is
    # unphased,
    print(bgen["X"][8, 3].compute().sel(data="phased").item())
    # the estimated probabilities imply the dosage
    # (or expected number of alleles)
    print(convert_to_dosage(p, nalleles, ploidy))
