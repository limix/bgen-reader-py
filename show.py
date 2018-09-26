from texttable import Texttable

from bgen_reader import (
    read_bgen,
    allele_expectation,
    example_files,
    compute_dosage,
    allele_frequency,
)

sampleid = "sample_005"
rsid = "RSID_6"

with example_files("example.32bits.bgen") as filepath:
    bgen = read_bgen(filepath, verbose=False)

    locus = bgen["variants"].query("rsid == '{}'".format(rsid)).index[0]
    sample = bgen["samples"].query("id == '{}'".format(sampleid)).index[0]

    nalleles = bgen["variants"].loc[locus, "nalleles"].item()
    ploidy = 2

    p = bgen["genotype"][locus, sample].compute()
    # For unphased genotypes only.
    e = allele_expectation(bgen["genotype"][locus, sample], nalleles, ploidy)

    alleles = bgen["variants"].loc[locus, "allele_ids"].split(",")

    tab = Texttable()

    tab.add_rows(
        [
            ["", "AA", "AG", "GG", "E[.]"],
            ["p"] + list(p) + [1.0],
            ["#" + alleles[0], 2, 1, 0, e[0]],
            ["#" + alleles[1], 0, 1, 2, e[1]],
        ]
    )

print(tab.draw())
print("variant: {}".format(rsid))
print("sample : {}".format(sampleid))
print()

e = allele_expectation(bgen["genotype"], nalleles, ploidy)

freq = allele_frequency(e)[locus]
print("Frequency of locus {}:".format(rsid))
print("    {}: {:f}".format(alleles[0], freq[0]))
print("    {}: {:f}".format(alleles[1], freq[1]))

# Alleles with minor allele frequencies accordong to the provided expections are used
# references by default.
dos = compute_dosage(e)
print()
print("Dosage: {:f}".format(dos[locus, sample]))
print()
