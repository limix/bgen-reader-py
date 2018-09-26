from bgen_reader import read_bgen, allele_expectation, example_files, compute_dosage

with example_files("example.32bits.bgen") as filepath:
    bgen = read_bgen(filepath, verbose=False)
    e = allele_expectation(bgen["genotype"], nalleles=2, ploidy=2)
    dos = compute_dosage(e)
    print(dos.shape)
    print(dos)
