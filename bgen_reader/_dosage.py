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
        Allele expectations encoded as a samples-by-alleles matrix.
    alt : array_like, optional
        Alternative allele index. If ``None``, the allele having the minor
        allele frequency for the provided ``expec`` is used as the alternative.
        Defaults to ``None``.

    Returns
    -------
    :class:`numpy.ndarray`
        Dosage encoded as an array of size equal to the number of samples.

    Examples
    --------
    .. doctest::

        >>> from pandas import DataFrame, option_context
        >>> from xarray import DataArray
        ...
        >>> from bgen_reader import (
        ...     allele_expectation,
        ...     allele_frequency,
        ...     compute_dosage,
        ...     example_files,
        ...     read_bgen,
        ... )
        ...
        >>> with example_files("example.32bits.bgen") as filepath: # doctest: +NORMALIZE_WHITESPACE +IGNORE_EXCEPTION_DETAIL
        ...     with option_context("display.max_rows", 6):
        ...         bgen = read_bgen(filepath, verbose=False)
        ...         variants = bgen["variants"]
        ...         genotype = bgen["genotype"]
        ...         samples = bgen["samples"]
        ...
        ...         variant_idx = 3
        ...         variant = variants.loc[variant_idx].compute()
        ...         print("---- Variant ----")
        ...         print(variant)
        ...         print()
        ...
        ...         geno = bgen["genotype"][variant_idx].compute()
        ...         metageno = DataFrame({k: geno[k] for k in ["ploidy", "missing"]},
        ...                              index=samples)
        ...         metageno.index.name = "sample"
        ...         print(metageno)
        ...         print()
        ...
        ...         p = DataArray(
        ...             geno["probs"],
        ...             name="probability",
        ...             coords={"sample": samples},
        ...             dims=["sample", "genotype"],
        ...         )
        ...         print("---- Probability ----")
        ...         print(p.to_series().unstack(level=-1))
        ...         print()
        ...
        ...         alleles = variant["allele_ids"].item().split(",")
        ...         e = DataArray(
        ...             allele_expectation(bgen, variant_idx),
        ...             name="expectation",
        ...             coords={"sample": samples, "allele": alleles},
        ...             dims=["sample", "allele"],
        ...         )
        ...         print("---- Expectation ----")
        ...         print(e.to_series().unstack(level=-1))
        ...         print()
        ...
        ...         rsid = variant["rsid"].item()
        ...         chrom = variant["chrom"].item()
        ...         variant_name = f"{chrom}:{rsid}"
        ...         f = DataFrame(allele_frequency(e), columns=[variant_name],
        ...                       index=alleles)
        ...         f.index.name = "allele"
        ...         print("---- Frequency ----")
        ...         print(f)
        ...         print()
        ...
        ...         alt = f.idxmin().item()
        ...         alt_idx = alleles.index(alt)
        ...         d = compute_dosage(e, alt=alt_idx).to_series()
        ...         d = DataFrame(d.values, columns=[f"alt={alt}"], index=d.index)
        ...         print("---- Dosage ----")
        ...         print(d)
        ...         print()
        ---- Variant ----
                id    rsid chrom   pos  nalleles allele_ids  vaddr
        3  SNPID_5  RSID_5    01  5000         2        A,G  16034
        <BLANKLINE>
                    ploidy  missing
        sample
        sample_001       2    False
        sample_002       2    False
        sample_003       2    False
        ...            ...      ...
        sample_498       2    False
        sample_499       2    False
        sample_500       2    False
        <BLANKLINE>
        [500 rows x 2 columns]
        <BLANKLINE>
        ---- Probability ----
        genotype           0         1         2
        sample
        sample_001  0.004883  0.028381  0.966736
        sample_002  0.990448  0.009277  0.000275
        sample_003  0.989319  0.003906  0.006775
        ...              ...       ...       ...
        sample_498  0.005524  0.994232  0.000244
        sample_499  0.012665  0.011536  0.975800
        sample_500  0.000214  0.984314  0.015472
        <BLANKLINE>
        [500 rows x 3 columns]
        <BLANKLINE>
        ---- Expectation ----
        allele             A         G
        sample
        sample_001  0.038147  1.961853
        sample_002  1.990173  0.009827
        sample_003  1.982544  0.017456
        ...              ...       ...
        sample_498  1.005280  0.994720
        sample_499  0.036865  1.963135
        sample_500  0.984742  1.015258
        <BLANKLINE>
        [500 rows x 2 columns]
        <BLANKLINE>
        ---- Frequency ----
                 01:RSID_5
        allele
        A       305.972181
        G       194.027819
        <BLANKLINE>
        ---- Dosage ----
                       alt=G
        sample
        sample_001  1.961853
        sample_002  0.009827
        sample_003  0.017456
        ...              ...
        sample_498  0.994720
        sample_499  1.963135
        sample_500  1.015258
        <BLANKLINE>
        [500 rows x 1 columns]
        <BLANKLINE>
    """
    if alt is None:
        return expec[..., -1]
    try:
        return expec[:, alt]
    except NotImplementedError:
        alt = asarray(alt, int)
        return asarray(expec, float)[:, alt]


def allele_expectation(bgen, variant_idx):
    r""" Allele expectation.

    Compute the expectation of each allele from the genotype probabilities.

    Parameters
    ----------
    bgen : bgen_file
        Bgen file handler.
    variant_idx : int
        Variant index.

    Returns
    -------
    :class:`numpy.ndarray`
        Samples-by-alleles matrix of allele expectations.

    Note
    ----
    This function supports unphased genotypes only.

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
