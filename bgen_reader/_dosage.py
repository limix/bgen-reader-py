from ._helper import get_genotypes, genotypes_to_allele_counts
from numpy import asarray


def convert_to_dosage(p, nalleles, ploidy):
    g = get_genotypes(ploidy, nalleles)
    c = genotypes_to_allele_counts(g)
    return (asarray(c, float).T * p).sum(1)
