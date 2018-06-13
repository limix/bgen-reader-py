from bgen_reader import read_bgen

if __name__ == "__main__":

    bgen = read_bgen("example.bgen", verbose=False)

    print(bgen["variants"].head())
    print(bgen["samples"].head())
    print(len(bgen["genotype"]))
    p = bgen["genotype"][0].compute()
    print(p)
    print(p.shape)
