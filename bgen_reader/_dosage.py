from numpy import asarray, stack

from ._helper import genotypes_to_allele_counts, get_genotypes


def allele_frequency(expec):
    r""" Compute allele frequency from its expectation.

    Parameters
    ----------
    expec : array_like
        Allele expectations encoded as a samples-by-alleles matrix.

    Returns
    -------
    :class:`numpy.ndarray`
        Allele frequencies encoded as a variants-by-alleles matrix.

    Examples
    --------
    .. doctest::

    >>> from texttable import Texttable
    >>> from bgen_reader import read_bgen, example_files
    >>> from bgen_reader import allele_expectation, allele_frequency
    ...
    >>> rsid = "RSID_6"
    ...
    >>> with example_files("example.32bits.bgen") as filepath:
    ...     bgen = read_bgen(filepath, verbose=False)
    ...     variants = bgen["variants"]
    ...     samples = bgen["samples"]
    ...     genotype = bgen["genotype"]
    ...
    ...     variant = variants[variants["rsid"] == rsid].compute()
    ...     variant_idx = variant.index.item()
    ...
    ...     p = genotype[variant_idx].compute()["probs"]
    ...     # For unphased genotypes only.
    ...     e = allele_expectation(bgen, variant_idx)
    ...     f = allele_frequency(e)
    ...
    ...     alleles = variant["allele_ids"].item().split(",")
    >>> print(alleles[0] + ": {}".format(f[0]))
    A: 229.23103218810434
    >>> print(alleles[1] + ": {}".format(f[1]))
    G: 270.7689678118956
    >>> print(variant)
            id    rsid chrom   pos  nalleles allele_ids  vaddr
    4  SNPID_6  RSID_6    01  6000         2        A,G  19377
    """
    expec = asarray(expec, float)
    if expec.ndim != 2:
        raise ValueError("Expectation matrix must be bi-dimensional.")
    ploidy = expec.shape[-1]
    return expec.sum(-2) / ploidy


def compute_dosage(expec, alt=None):
    r""" Compute dosage from allele expectation.

    Parameters
    ----------
    expec : array_like
        Allele expectations encoded as a variants-by-samples-by-alleles matrix.
    ref : array_like
        Allele reference of each locus. The allele having the minor allele frequency
        for the provided ``expec`` is used as the reference if `None`. Defaults to
        `None`.

    Returns
    -------
    :class:`numpy.ndarray`
        Dosage encoded as a variants-by-samples matrix.

    Examples
    --------

    .. doctest::
    >>> from bgen_reader import read_bgen, allele_expectation, example_files
    >>> from bgen_reader import compute_dosage
    ...
    >>> with example_files("example.32bits.bgen") as filepath:
    ...     bgen = read_bgen(filepath, verbose=False)
    ...     e = allele_expectation(bgen, 0)
    ...     dosage = compute_dosage(e)
    ...     print(dosage.shape)
    ...     print(dosage[:5])
    (500,)
    [       nan 1.93575854 1.91558579 1.0174256  1.91159064]
    """
    if alt is None:
        return expec[..., -1]
    try:
        return expec[alt, :, alt]
    except NotImplementedError:
        alt = asarray(alt, int)
        return asarray(expec, float)[alt, :, alt]


def allele_expectation(bgen, variant_idx):
    r""" Allele expectation.

    Compute the expectation of each allele from the bgen probabilities.

    Parameters
    ----------
    bgen : bgen_file
        Bgen file handler.
    variant_idx : int
        Variant index.

    Returns
    -------
    :class:`numpy.ndarray`
        Samples by allele matrix of ellele expectation.

    Examples
    --------
    .. doctest::

    >>> from texttable import Texttable
    >>> from bgen_reader import read_bgen, allele_expectation, example_files
    ...
    >>> sampleid = "sample_005"
    >>> rsid = "RSID_6"
    ...
    >>> with example_files("example.32bits.bgen") as filepath:
    ...     bgen = read_bgen(filepath, verbose=False)
    ...     variants = bgen["variants"]
    ...     samples = bgen["samples"]
    ...     genotype = bgen["genotype"]
    ...
    ...     variant = variants[variants["rsid"] == rsid].compute()
    ...     variant_idx = variant.index.item()
    ...
    ...     sample_idx = samples[samples == sampleid].index.item()
    ...
    ...     p = genotype[variant_idx].compute()["probs"][sample_idx]
    ...     # For unphased genotypes only.
    ...     e = allele_expectation(bgen, variant_idx)[sample_idx]
    ...
    ...     alleles = variant["allele_ids"].item().split(",")
    ...
    ...     print(Texttable().add_rows(
    ...         [
    ...             ["", "AA", "AG", "GG", "E[.]"],
    ...             ["p"] + list(p) + ["na"],
    ...             ["#" + alleles[0], 2, 1, 0, e[0]],
    ...             ["#" + alleles[1], 0, 1, 2, e[1]],
    ...         ]
    ...     ).draw())
    +----+-------+-------+-------+-------+
    |    |  AA   |  AG   |  GG   | E[.]  |
    +====+=======+=======+=======+=======+
    | p  | 0.012 | 0.987 | 0.001 | na    |
    +----+-------+-------+-------+-------+
    | #A | 2     | 1     | 0     | 1.011 |
    +----+-------+-------+-------+-------+
    | #G | 0     | 1     | 2     | 0.989 |
    +----+-------+-------+-------+-------+
    >>> print(e)
    [1.01086423 0.98913577]
    >>> print(variant)
            id    rsid chrom   pos  nalleles allele_ids  vaddr
    4  SNPID_6  RSID_6    01  6000         2        A,G  19377

    Note
    ----
    This function supports unphased genotypes only.
    """
    geno = bgen["genotype"][variant_idx].compute()
    if geno["phased"]:
        raise ValueError("Allele expectation is define for unphased genotypes only.")

    nalleles = bgen["variants"].loc[variant_idx, "nalleles"].compute().item()
    genotypes = get_genotypes(geno["ploidy"], nalleles)
    expec = []
    for i in range(len(genotypes)):
        count = asarray(genotypes_to_allele_counts(genotypes[i]), float)
        n = count.shape[0]
        expec.append((count.T * geno["probs"][i, :n]).sum(1))

    return stack(expec, axis=0)
