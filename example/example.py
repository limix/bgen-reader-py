from bgen_reader import read_bgen

(variants, samples, genotype) = read_bgen("example.bgen")

print(variants.head())
print(samples.head())
print(len(genotype))
print(genotype[0].compute())
