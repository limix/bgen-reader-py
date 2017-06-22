from bgen_reader import read_bgen

bgen = read_bgen("example.bgen", verbose=False)

print(bgen['variants'].head())
print(bgen['samples'].head())
print(len(bgen['genotype']))
print(bgen['genotype'][0].compute())
