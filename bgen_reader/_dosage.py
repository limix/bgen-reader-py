from ._helper import get_genotypes, genotypes_to_allele_counts
from numpy import asarray, newaxis


def convert_to_dosage(p, nalleles, ploidy):
    r"""Convert probabilities to dosage.

    Parameters
    ----------
    p : array_like
        Allele probabilities.
    nalleles : int
        Number of alleles.
    ploidy : int
        Number of complete sets of chromosomes.

    Returns
    -------
    :class:`numpy.ndarray`
        Dosage matrix.

    Warning
    -------
    This is a new function that needs more testing. Please, report any problem.
    """
    g = get_genotypes(ploidy, nalleles)
    c = genotypes_to_allele_counts(g)
    return (asarray(c, float).T * p).sum(1)


def allele_expectation(p, nalleles, ploidy):
    r"""Allele expectation.

    Compute the expectation of each allele from the given probabilities.

    Parameters
    ----------
    p : array_like
        Allele probabilities.
    nalleles : int
        Number of alleles.
    ploidy : int
        Number of complete sets of chromosomes.

    Returns
    -------
    :class:`numpy.ndarray`
        Expectation of each allele.

    Warning
    -------
    This is a new function that needs more testing. Please, report any problem.
    """
    import pdb

    g = get_genotypes(ploidy, nalleles)
    c = genotypes_to_allele_counts(g)
    c = asarray(c, float)

    # return (asarray(c, float).T * p.compute()).sum(1)
    # return (p.compute() * c.reshape((c.shape[0], 1, 1, c.shape[1]))).sum(1)
    r = c.reshape((1,) * (p.ndim - 1) + (c.shape[1], c.shape[0]))
    return (r * p.compute()[..., newaxis, :]).sum(-1)
