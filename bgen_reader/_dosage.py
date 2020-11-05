from numpy import asarray, stack

from ._helper import genotypes_to_allele_counts, get_genotypes


def allele_frequency(expec):
    r"""Compute allele frequency from its expectation.

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

        >>> from bgen_reader import read_bgen, example_filepath
        >>> from bgen_reader import allele_expectation, allele_frequency
        >>>
        >>> filepath = example_filepath("example.32bits.bgen")
        >>>
        >>> bgen = read_bgen(filepath, verbose=False)
        >>>
        >>> variants = bgen["variants"]
        >>> samples = bgen["samples"]
        >>> genotype = bgen["genotype"]
        >>>
        >>> variant = variants[variants["rsid"] == "RSID_6"].compute()
        >>> variant_idx = variant.index.values[0]
        >>>
        >>> p = genotype[variant_idx].compute()["probs"]
        >>> # For unphased genotypes only.
        >>> e = allele_expectation(bgen, variant_idx)
        >>> f = allele_frequency(e)
        >>>
        >>> alleles = variant["allele_ids"].values[0].split(",")
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
    nallele0 = expec.shape[-1]
    return expec.sum(-2) / nallele0


def compute_dosage(expec, alt=None):
    r"""Compute dosage from allele expectation.

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
    .. code-block:: python
        :caption: First a quick-start example.

        >>> from bgen_reader import allele_expectation, compute_dosage
        >>> from bgen_reader import example_filepath, read_bgen
        >>>
        >>> filepath = example_filepath("example.32bits.bgen")
        >>>
        >>> # Read the example.
        >>> bgen = read_bgen(filepath, verbose=False)
        >>>
        >>> # Extract the allele expectations of the fourth variant.
        >>> variant_idx = 3
        >>> e = allele_expectation(bgen, variant_idx)
        >>>
        >>> # Compute the dosage when considering the first allele
        >>> # as the reference/alternative one.
        >>> alt_allele_idx = 1
        >>> d = compute_dosage(e, alt=alt_allele_idx)
        >>>
        >>> # Print the dosage of the first five samples only.
        >>> print(d[:5])
        [1.96185308 0.00982666 0.01745552 1.00347899 1.01153563]

    .. code-block:: python
        :caption: Genotype probabilities, allele expectations and frequencies.

        >>> from bgen_reader import (
        ...     allele_expectation,
        ...     allele_frequency,
        ...     compute_dosage,
        ...     example_filepath,
        ...     read_bgen,
        ... )
        >>> from pandas import DataFrame
        >>> from xarray import DataArray
        >>>
        >>> # Download an example
        >>> filepath = example_filepath("example.32bits.bgen")
        >>>
        >>> # Open the bgen file.
        >>> bgen = read_bgen(filepath, verbose=False)
        >>> variants = bgen["variants"]
        >>> genotype = bgen["genotype"]
        >>> samples = bgen["samples"]
        >>>
        >>> variant_idx = 3
        >>> variant = variants.loc[variant_idx].compute()
        >>> # Print the metadata of the fourth variant.
        >>> print(variant)
                id    rsid chrom   pos  nalleles allele_ids  vaddr
        3  SNPID_5  RSID_5    01  5000         2        A,G  16034

        >>> geno = bgen["genotype"][variant_idx].compute()
        >>> metageno = DataFrame({k: geno[k] for k in ["ploidy", "missing"]},
        ...                      index=samples)
        >>> metageno.index.name = "sample"
        >>> print(metageno)  # doctest: +NORMALIZE_WHITESPACE
                    ploidy  missing
        sample
        sample_001       2    False
        sample_002       2    False
        sample_003       2    False
        sample_004       2    False
        ...            ...      ...
        sample_497       2    False
        sample_498       2    False
        sample_499       2    False
        sample_500       2    False
        <BLANKLINE>
        [500 rows x 2 columns]
        >>> p = DataArray(
        ...     geno["probs"],
        ...     name="probability",
        ...     coords={"sample": samples},
        ...     dims=["sample", "genotype"],
        ... )
        >>> # Print the genotype probabilities.
        >>> print(p.to_series().unstack(level=-1))  # doctest: +NORMALIZE_WHITESPACE
        genotype          0        1        2
        sample
        sample_001  0.00488  0.02838  0.96674
        sample_002  0.99045  0.00928  0.00027
        sample_003  0.98932  0.00391  0.00677
        sample_004  0.00662  0.98328  0.01010
        ...             ...      ...      ...
        sample_497  0.00137  0.01312  0.98550
        sample_498  0.00552  0.99423  0.00024
        sample_499  0.01266  0.01154  0.97580
        sample_500  0.00021  0.98431  0.01547
        <BLANKLINE>
        [500 rows x 3 columns]
        >>> alleles = variant["allele_ids"].values[0].split(",")
        >>> e = DataArray(
        ...     allele_expectation(bgen, variant_idx),
        ...     name="expectation",
        ...     coords={"sample": samples, "allele": alleles},
        ...     dims=["sample", "allele"],
        ... )
        >>> # Print the allele expectations.
        >>> print(e.to_series().unstack(level=-1))  # doctest: +NORMALIZE_WHITESPACE
        allele            A        G
        sample
        sample_001  0.03815  1.96185
        sample_002  1.99017  0.00983
        sample_003  1.98254  0.01746
        sample_004  0.99652  1.00348
        ...             ...      ...
        sample_497  0.01587  1.98413
        sample_498  1.00528  0.99472
        sample_499  0.03687  1.96313
        sample_500  0.98474  1.01526
        <BLANKLINE>
        [500 rows x 2 columns]
        >>> rsid = variant["rsid"].values[0]
        >>> chrom = variant["chrom"].values[0]
        >>> variant_name = f"{chrom}:{rsid}"
        >>> f = DataFrame(allele_frequency(e), columns=[variant_name], index=alleles)
        >>> f.index.name = "allele"
        >>> # Allele frequencies.
        >>> print(f) # doctest: +NORMALIZE_WHITESPACE
                01:RSID_5
        allele
        A       305.97218
        G       194.02782
        >>> alt = f.idxmin().values[0]
        >>> alt_idx = alleles.index(alt)
        >>> d = compute_dosage(e, alt=alt_idx).to_series()
        >>> d = DataFrame(d.values, columns=[f"alt={alt}"], index=d.index)
        >>> # Dosages when considering G as the alternative allele.
        >>> print(d) # doctest: +NORMALIZE_WHITESPACE
                      alt=G
        sample
        sample_001  1.96185
        sample_002  0.00983
        sample_003  0.01746
        sample_004  1.00348
        ...             ...
        sample_497  1.98413
        sample_498  0.99472
        sample_499  1.96313
        sample_500  1.01526
        <BLANKLINE>
        [500 rows x 1 columns]
    """
    if alt is None:
        return expec[..., -1]
    try:
        return expec[:, alt]
    except NotImplementedError:
        alt = asarray(alt, int)
        return asarray(expec, float)[:, alt]


def allele_expectation(bgen, variant_idx):
    r"""Allele expectation.

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

        >>> from bgen_reader import allele_expectation, example_filepath, read_bgen
        >>>
        >>> from texttable import Texttable
        >>>
        >>> filepath = example_filepath("example.32bits.bgen")
        >>>
        >>> # Read the example.
        >>> bgen = read_bgen(filepath, verbose=False)
        >>>
        >>> variants = bgen["variants"]
        >>> samples = bgen["samples"]
        >>> genotype = bgen["genotype"]
        >>>
        >>> genotype = bgen["genotype"]
        >>> # This `compute` call will return a pandas data frame,
        >>> variant = variants[variants["rsid"] == "RSID_6"].compute()
        >>> # from which we retrieve the variant index.
        >>> variant_idx = variant.index.values[0]
        >>> print(variant)
                id    rsid chrom   pos  nalleles allele_ids  vaddr
        4  SNPID_6  RSID_6    01  6000         2        A,G  19377
        >>> genotype = bgen["genotype"]
        >>> # Samples is a pandas series, and we retrieve the
        >>> # sample index from the sample name.
        >>> sample_idx = samples[samples == "sample_005"].index.values[0]
        >>>
        >>> genotype = bgen["genotype"]
        >>> # This `compute` call will return a dictionary from which
        >>> # we can get the probability matrix the corresponding
        >>> # variant.
        >>> p = genotype[variant_idx].compute()["probs"][sample_idx]
        >>>
        >>> genotype = bgen["genotype"]
        >>> # Allele expectation makes sense for unphased genotypes only,
        >>> # which is the case here.
        >>> e = allele_expectation(bgen, variant_idx)[sample_idx]
        >>>
        >>> genotype = bgen["genotype"]
        >>> alleles = variant["allele_ids"].values[0].split(",")
        >>>
        >>> genotype = bgen["genotype"]
        >>>
        >>> # Print what we have got in a nice format.
        >>> table = Texttable()
        >>> table = table.add_rows(
        ...     [
        ...         ["", "AA", "AG", "GG", "E[.]"],
        ...         ["p"] + list(p) + ["na"],
        ...         ["#" + alleles[0], 2, 1, 0, e[0]],
        ...         ["#" + alleles[1], 0, 1, 2, e[1]],
        ...     ]
        ... )
        >>> print(table.draw())
        +----+-------+-------+-------+-------+
        |    |  AA   |  AG   |  GG   | E[.]  |
        +====+=======+=======+=======+=======+
        | p  | 0.012 | 0.987 | 0.001 | na    |
        +----+-------+-------+-------+-------+
        | #A | 2     | 1     | 0     | 1.011 |
        +----+-------+-------+-------+-------+
        | #G | 0     | 1     | 2     | 0.989 |
        +----+-------+-------+-------+-------+
    """
    geno = bgen["genotype"][variant_idx].compute()
    if geno["phased"]:
        raise ValueError("Allele expectation is define for unphased genotypes only.")

    nalleles = bgen["variants"].loc[variant_idx, "nalleles"].compute().values[0]
    genotypes = get_genotypes(geno["ploidy"], nalleles)
    expec = []
    for i, genotype in enumerate(genotypes):
        count = asarray(genotypes_to_allele_counts(genotype), float)
        n = count.shape[0]
        expec.append((count.T * geno["probs"][i, :n]).sum(1))

    return stack(expec, axis=0)
