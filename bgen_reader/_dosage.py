from ._helper import get_genotypes, genotypes_to_allele_counts
from numpy import asarray, newaxis, where, ones, empty, stack, arange, ix_


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


def allele_frequency(expec):
    ploidy = expec.shape[-1]
    if expec.ndim < 3:
        n = 1
    else:
        n = expec.shape[1]
    return expec.sum(-2) / (ploidy * n)


def compute_dosage(expec, ref=None):
    freq = allele_frequency(expec)
    if ref is None:
        ma = freq.argmin(1)
    ma = asarray(ma, int)
    return expec[ma, :, ma]


def allele_expectation(p, nalleles, ploidy):
    r"""Allele expectation.

    Compute the expectation of each allele from the given probabilities.
    It accepts three shapes of matrices:
    - unidimensional probability array;
    - bidimensional sample-by-probability array; and
    - three dimensional locus-by-sample-by-probability array.

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
        Last dimension will contain the expectation of each allele.

    Examples
    --------

    .. codetest::

    >>> from texttable import Texttable
    >>> from bgen_reader import read_bgen, allele_expectation, example_files
    >>>
    >>> sampleid = "sample_005"
    >>> rsid = "RSID_6"
    >>>
    >>> with example_files("example.32bits.bgen") as filepath:
    ...     bgen = read_bgen(filepath, verbose=False)
    ...
    ...     locus = bgen["variants"].query("rsid == '{}'".format(rsid)).index
    ...     sample = bgen["samples"].query("id == '{}'".format(sampleid)).index
    ...
    ...     nalleles = bgen["variants"].loc[locus, "nalleles"].item()
    ...     ploidy = 2
    ...
    ...     p = bgen["genotype"][locus[0], sample[0]].compute()
    ...     # For unphased genotypes only.
    ...     e = allele_expectation(bgen["genotype"][locus[0], sample[0]], nalleles, ploidy)
    ...
    ...     alleles = bgen["variants"].loc[locus, "allele_ids"].item().split(",")
    ...
    ...     tab = Texttable()
    ...
    ...     tab.add_rows(
    ...         [
    ...             ["", "AA", "AG", "GG", "E[.]"],
    ...             ["p"] + list(p) + [1.0],
    ...             ["#" + alleles[0], 2, 1, 0, e[0]],
    ...             ["#" + alleles[1], 0, 1, 2, e[1]],
    ...         ]
    ...     )
    >>> print(tab.draw())
    +----+-------+-------+-------+-------+
    |    |  AA   |  AG   |  GG   | E[.]  |
    +====+=======+=======+=======+=======+
    | p  | 0.012 | 0.987 | 0.001 | 1     |
    +----+-------+-------+-------+-------+
    | #A | 2     | 1     | 0     | 1.011 |
    +----+-------+-------+-------+-------+
    | #G | 0     | 1     | 2     | 0.989 |
    +----+-------+-------+-------+-------+
    >>> print("variant: {}".format(rsid))
    variant: RSID_6
    >>> print("sample : {}".format(sampleid))
    sample : sample_005

    Warning
    -------
    This function supports unphased genotypes only.
    """
    g = get_genotypes(ploidy, nalleles)
    c = asarray(genotypes_to_allele_counts(g), float)
    c = c.T.reshape((1,) * (p.ndim - 1) + (c.shape[1], c.shape[0]))
    return (c * p.compute()[..., newaxis, :]).sum(-1)
