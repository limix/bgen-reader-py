***********
Quick-start
***********

The following examples assume you have downloaded the ``example.bgen``,
``haplotypes.bgen``, and ``complex.bgen`` files (found in this repository) to the
directory you are executing Python.

.. doctest::

    >>> from bgen_reader import download
    >>>
    >>> download("http://rest.s3for.me/bgen-reader/complex.bgen", verbose=False)
    >>> download("http://rest.s3for.me/bgen-reader/example.bgen", verbose=False)
    >>> download("http://rest.s3for.me/bgen-reader/haplotypes.bgen", verbose=False)

Unphased genotype
=================

.. doctest::

   >>> from bgen_reader import read_bgen
   >>>
   >>> bgen = read_bgen("example.bgen", verbose=False)
   >>>
   >>> print(bgen["variants"].head())
           id    rsid chrom   pos  nalleles allele_ids
   0  SNPID_2  RSID_2    01  2000         2        A,G
   1  SNPID_3  RSID_3    01  3000         2        A,G
   2  SNPID_4  RSID_4    01  4000         2        A,G
   3  SNPID_5  RSID_5    01  5000         2        A,G
   4  SNPID_6  RSID_6    01  6000         2        A,G

.. doctest::

   >>> print(bgen["samples"].head())
              id
   0  sample_001
   1  sample_002
   2  sample_003
   3  sample_004
   4  sample_005

.. doctest::

   >>> print(len(bgen["genotype"]))
   199

.. doctest::

   >>> p = bgen["genotype"][0].compute()
   >>> print(p)
   [[       nan        nan        nan]
    [0.02780236 0.00863674 0.9635609 ]
    [0.01736504 0.04968414 0.93295083]
    ...
    [0.01419069 0.02810669 0.95770262]
    [0.91949463 0.05206298 0.02844239]
    [0.00244141 0.98410029 0.0134583 ]]

.. doctest::

   >>> print(p.shape)
   (500, 3)

The ``example.bgen`` file can be found in the ``example`` folder, as
well as the next ones.


Phased genotype
===============

.. doctest::

   >>> from bgen_reader import read_bgen
   >>> bgen = read_bgen("haplotypes.bgen", verbose=False)
   >>>
   >>> print(bgen["variants"].head())
        id rsid chrom  pos  nalleles allele_ids
   0  SNP1  RS1     1    1         2        A,G
   1  SNP2  RS2     1    2         2        A,G
   2  SNP3  RS3     1    3         2        A,G
   3  SNP4  RS4     1    4         2        A,G

.. doctest::

   >>> print(bgen["samples"].head())
            id
   0  sample_0
   1  sample_1
   2  sample_2
   3  sample_3

.. doctest::

   >>> # Print the estimated probabilities for the first variant
   >>> # and second individual.
   >>> print(bgen["genotype"][0, 1].compute())
   [0. 1. 1. 0.]

.. doctest::

   >>> # Is it a phased one?
   >>> print(bgen["X"][0, 1].compute().sel(data="phased").item())
   1

.. doctest::

   >>> # How many haplotypes?
   >>> print(bgen["X"][0, 1].compute().sel(data="ploidy").item())
   2

.. doctest::

   >>> # And how many alleles?
   >>> print(bgen["variants"].loc[0, "nalleles"])
   2

.. doctest::

   >>> # Therefore, the first haplotype has probability 100%
   >>> # of having the allele
   >>> print(bgen["variants"].loc[0, "allele_ids"].split(",")[1])
   G

.. doctest::

   >>> # And the second haplotype has probability 100% of having
   >>> # the first allele
   >>> print(bgen["variants"].loc[0, "allele_ids"].split(",")[0])
   A

Complex file
============

.. doctest::

   >>> from bgen_reader import read_bgen
   >>>
   >>> bgen = read_bgen("complex.bgen", verbose=False)
   >>>
   >>> print(bgen["variants"])
        id rsid chrom  pos  nalleles                            allele_ids
   0         V1    01    1         2                                   A,G
   1  V2.1   V2    01    2         2                                   A,G
   2         V3    01    3         2                                   A,G
   3         M4    01    4         3                                 A,G,T
   4         M5    01    5         2                                   A,G
   5         M6    01    7         4                            A,G,GT,GTT
   6         M7    01    7         6                 A,G,GT,GTT,GTTT,GTTTT
   7         M8    01    8         7          A,G,GT,GTT,GTTT,GTTTT,GTTTTT
   8         M9    01    9         8  A,G,GT,GTT,GTTT,GTTTT,GTTTTT,GTTTTTT
   9        M10    01   10         2                                   A,G

.. doctest::

   >>> print(bgen["samples"])
            id
   0  sample_0
   1  sample_1
   2  sample_2
   3  sample_3

.. doctest::

   >>> # Print the estimated probabilities for the first variant
   >>> # and second individual.
   >>> print(bgen["genotype"][0, 1].compute())
   [ 1.  0.  0. nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan
    nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan]

.. doctest::

   >>> # The NaN elements are a by-product of the heterogenous
   >>> # ploidy and number of alleles across variants and samples.
   >>> # For example, the 9th variant for the 4th individual
   >>> # has ploidy
   >>> ploidy = bgen["X"][8, 3].compute().sel(data="ploidy").item()
   >>> print(ploidy)
   2

.. doctest::

   >>> # and number of alleles equal to
   >>> nalleles = bgen["variants"].loc[8, "nalleles"]
   >>> print(nalleles)
   8

.. doctest::

   >>> # Its probability distribution is given by the array
   >>> p = bgen["genotype"][8, 3].compute()
   >>> print(p)
   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
    0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]

.. doctest::

   >>> # of size
   >>> print(len(p))
   36

.. doctest::

   >>> # Since the 9th variant for the 4th individual is
   >>> # unphased,
   >>> print(bgen["X"][8, 3].compute().sel(data="phased").item())
   0

Dosage
======

For a genotype with ploidy two and locus with two possible alleles, the dosage
is defined as the expectation of the number of the reference alleles.
It is common to define the reference allele as being the one has lower frequency
under the given dataset.
The following example demonstrate that case.

.. doctest::

    >>> from bgen_reader import read_bgen, allele_expectation, example_files
    >>> from bgen_reader import compute_dosage
    >>>
    >>> with example_files("example.32bits.bgen") as filepath:
    ...     bgen = read_bgen(filepath, verbose=False)
    ...     e = allele_expectation(bgen["genotype"], nalleles=2, ploidy=2)
    ...     dosage = compute_dosage(e)
    ...     print(dosage.shape)
    ...     print(dosage)
    (199, 500)
    [[       nan 0.06424146 0.08441421 ... 0.05648808 1.89105224 0.98898311]
    [1.98779296 1.97802735 0.02111815 ... 1.95492412 1.00897216 1.02255316]
    [       nan 0.06424146 0.08441421 ... 0.05648808 1.89105224 0.98898311]
    ...
    [       nan 0.06424146 0.08441421 ... 0.05648808 1.89105224 0.98898311]
    [1.98779296 1.97802735 0.02111815 ... 1.95492412 1.00897216 1.02255316]
    [1.98779296 1.97802735 0.02111815 ... 1.95492412 1.00897216 1.02255316]]

The function `compute_dosage` also accepts the argument `ref` from which the reference
alleles can be specified. (Consult `help(bgen_reader.compute_dosage)` for the full
specification.)

Another example now querying specific locus and sample.

.. doctest::

   >>> from texttable import Texttable
   >>>
   >>> from bgen_reader import (
   >>>     read_bgen,
   >>>     allele_expectation,
   >>>     example_files,
   >>>     compute_dosage,
   >>>     allele_frequency,
   >>> )
   >>>
   >>> sampleid = "sample_005"
   >>> rsid = "RSID_6"
   >>>
   >>> with example_files("example.32bits.bgen") as filepath:
   ...     bgen = read_bgen(filepath, verbose=False)
   ...
   ...     locus = bgen["variants"].query("rsid == '{}'".format(rsid)).index[0]
   ...     sample = bgen["samples"].query("id == '{}'".format(sampleid)).index[0]
   ...
   ...     nalleles = bgen["variants"].loc[locus, "nalleles"].item()
   ...     ploidy = 2
   ...
   ...     p = bgen["genotype"][locus, sample].compute()
   ...     # For unphased genotypes only.
   ...     e = allele_expectation(bgen["genotype"][locus, sample], nalleles, ploidy)
   ...
   ...     alleles = bgen["variants"].loc[locus, "allele_ids"].split(",")
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
   >>>
   >>> print(tab.draw())
   >>> print("variant: {}".format(rsid))
   >>> print("sample : {}".format(sampleid))
   >>>
   >>> e = allele_expectation(bgen["genotype"], nalleles, ploidy)
   >>>
   >>> freq = allele_frequency(e)[locus]
   >>> print("Frequency of locus {}:".format(rsid))
   >>> print("    {}: {:f}".format(alleles[0], freq[0]))
   >>> print("    {}: {:f}".format(alleles[1], freq[1]))
   >>>
   >>> # Alleles with minor allele frequencies accordong to the provided expections are used
   >>> # references by default.
   >>> dos = compute_dosage(e)
   >>> print()
   >>> print("Dosage: {:f}".format(dos[locus, sample]))
   >>> print()
   +----+-------+-------+-------+-------+
   |    |  AA   |  AG   |  GG   | E[.]  |
   +====+=======+=======+=======+=======+
   | p  | 0.012 | 0.987 | 0.001 | 1     |
   +----+-------+-------+-------+-------+
   | #A | 2     | 1     | 0     | 1.011 |
   +----+-------+-------+-------+-------+
   | #G | 0     | 1     | 2     | 0.989 |
   +----+-------+-------+-------+-------+
   variant: RSID_6
   sample : sample_005

   Frequency of locus RSID_6:
       A: 0.458462
       G: 0.541538

   Dosage: 0.088409
