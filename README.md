# bgen-reader

[![Travis](https://img.shields.io/travis/limix/bgen-reader-py.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.org/limix/bgen-reader-py) [![AppVeyor](https://img.shields.io/appveyor/ci/Horta/bgen-reader-py.svg?style=flat-square&label=windows%20build)](https://ci.appveyor.com/project/Horta/bgen-reader-py)

A [BGEN file format](http://www.well.ox.ac.uk/~gav/bgen_format/) reader.

BGEN is a file format for storing large genetic datasets.
It supports both unphased genotypes and phased haplotype data with variable
ploidy and number of alleles. It was designed to provides a compact data
representation without sacrificing variant access performance.

This Python package is a wrapper around the [bgen library](https://github.com/limix/bgen),
a low-memory footprint reader that efficiently reads BGEN files.
It fully supports the BGEN format specifications: 1.2 and 1.3;
as well as their optional compressed formats.

## Table of Contents

* [Install](#install)
* [Usage](#usage)
  * [Unphased genotype](#unphased-genotype)
  * [Phased genotype](#phased-genotype)
  * [Complex file](#complex-file)
  * [Dosage](#dosage)
* [Troubleshooting](#troubleshooting)
  * [fatal error: bgen.h: No such file or directory](#fatal-error-bgenh-no-such-file-or-directory)
* [Problems](#problems)
* [Authors](#authors)
* [License](#license)

## Install

The recommended way to install this package is via [conda](https://conda.io/docs/)

```bash
conda install -c conda-forge bgen-reader
```

Alternatively, it can be installed using the pip command

```bash
pip install bgen-reader
```

However, this method will require that
the [bgen](https://github.com/limix/bgen) C library has
been installed before.

## Usage

The following examples assume you have downloaded the ``example.bgen``,
``haplotypes.bgen``, and ``complex.bgen`` files (found in this repository) to the
directory you are executing Python.

### Unphased genotype

```python
>>> from bgen_reader import read_bgen
>>>
>>> bgen = read_bgen("example.bgen", verbose=False)
>>>
>>> print(bgen["variants"].head())
```

```bash
        id    rsid chrom   pos  nalleles allele_ids
0  SNPID_2  RSID_2    01  2000         2        A,G
1  SNPID_3  RSID_3    01  3000         2        A,G
2  SNPID_4  RSID_4    01  4000         2        A,G
3  SNPID_5  RSID_5    01  5000         2        A,G
4  SNPID_6  RSID_6    01  6000         2        A,G
```

```python
>>> print(bgen["samples"].head())
```

```bash
           id
0  sample_001
1  sample_002
2  sample_003
3  sample_004
4  sample_005
```

```python
>>> print(len(bgen["genotype"]))
```

```bash
199
```

```python
>>> p = bgen["genotype"][0].compute()
>>> print(p)
```

```bash
[[       nan        nan        nan]
 [0.02780236 0.00863674 0.9635609 ]
 [0.01736504 0.04968414 0.93295083]
 ...
 [0.01419069 0.02810669 0.95770262]
 [0.91949463 0.05206298 0.02844239]
 [0.00244141 0.98410029 0.0134583 ]]
```

```python
>>> print(p.shape)
```

```bash
(500, 3)
```

The ``example.bgen`` file can be found in the ``example`` folder, as
well as the next ones.

### Phased genotype

```python
>>> from bgen_reader import read_bgen
>>> bgen = read_bgen("haplotypes.bgen", verbose=False)
>>>
>>> print(bgen["variants"].head())
```

```bash
     id rsid chrom  pos  nalleles allele_ids
0  SNP1  RS1     1    1         2        A,G
1  SNP2  RS2     1    2         2        A,G
2  SNP3  RS3     1    3         2        A,G
3  SNP4  RS4     1    4         2        A,G
```

```python
>>> print(bgen["samples"].head())
```

```bash
         id
0  sample_0
1  sample_1
2  sample_2
3  sample_3
```

```python
>>> # Print the estimated probabilities for the first variant
>>> # and second individual.
>>> print(bgen["genotype"][0, 1].compute())
```
```bash
[0. 1. 1. 0.]
```

```python
>>> # Is it a phased one?
>>> print(bgen["X"][0, 1].compute().sel(data="phased").item())
```

```bash
1
```

```python
>>> # How many haplotypes?
>>> print(bgen["X"][0, 1].compute().sel(data="ploidy").item())
```

```bash
2
```

```python
>>> # And how many alleles?
>>> print(bgen["variants"].loc[0, "nalleles"])
```

```bash
2
```

```python
>>> # Therefore, the first haplotype has probability 100%
>>> # of having the allele
>>> print(bgen["variants"].loc[0, "allele_ids"].split(",")[1])
```

```bash
G
```

```python
>>> # And the second haplotype has probability 100% of having
>>> # the first allele
>>> print(bgen["variants"].loc[0, "allele_ids"].split(",")[0])
```

```bash
A
```

### Complex file

```python
>>> from bgen_reader import read_bgen, convert_to_dosage
>>>
>>> bgen = read_bgen("complex.bgen", verbose=False)
>>>
>>> print(bgen["variants"])
```

```bash
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
```

```python
>>> print(bgen["samples"])
```

```bash
         id
0  sample_0
1  sample_1
2  sample_2
3  sample_3
```

```python
>>> # Print the estimated probabilities for the first variant
>>> # and second individual.
>>> print(bgen["genotype"][0, 1].compute())
```

```bash
[ 1.  0.  0. nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan
 nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan]
```

```python
>>> # The NaN elements are a by-product of the heterogenous
>>> # ploidy and number of alleles across variants and samples.
>>> # For example, the 9th variant for the 4th individual
>>> # has ploidy
>>> ploidy = bgen["X"][8, 3].compute().sel(data="ploidy").item()
>>> print(ploidy)
```

```bash
2
```

```python
>>> # and number of alleles equal to
>>> nalleles = bgen["variants"].loc[8, "nalleles"]
>>> print(nalleles)
```

```
8
```

```python
>>> # Its probability distribution is given by the array
>>> p = bgen["genotype"][8, 3].compute()
>>> print(p)
```

```bash
[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
```

```python
>>> # of size
>>> print(len(p))
```

```bash
36
```

```python
>>> # Since the 9th variant for the 4th individual is
>>> # unphased,
>>> print(bgen["X"][8, 3].compute().sel(data="phased").item())
```

```bash
0
```

```python
>>> # the estimated probabilities imply the dosage
>>> # (or expected number of alleles)
>>> print(convert_to_dosage(p, nalleles, ploidy))
```

```bash
[0. 1. 0. 0. 0. 1. 0. 0.]
```

### Dosage

For a genotype with ploidy two and locus with two possible alleles, the dosage
is defined as the expectation of the number of the reference alleles.
It is common to define the reference allele as being the one has lower frequency
under the given dataset.
The following example demonstrate that case.

```python
>>> from bgen_reader import read_bgen, allele_expectation, example_files
>>> from bgen_reader import compute_dosage
>>>
>>> with example_files("example.32bits.bgen") as filepath:
...     bgen = read_bgen(filepath, verbose=False)
...     e = allele_expectation(bgen["genotype"], nalleles=2, ploidy=2)
...     dosage = compute_dosage(e)
...     print(dosage.shape)
...     print(dosage)
```

```bash
(199, 500)
[[       nan 0.06424146 0.08441421 ... 0.05648808 1.89105224 0.98898311]
[1.98779296 1.97802735 0.02111815 ... 1.95492412 1.00897216 1.02255316]
[       nan 0.06424146 0.08441421 ... 0.05648808 1.89105224 0.98898311]
...
[       nan 0.06424146 0.08441421 ... 0.05648808 1.89105224 0.98898311]
[1.98779296 1.97802735 0.02111815 ... 1.95492412 1.00897216 1.02255316]
[1.98779296 1.97802735 0.02111815 ... 1.95492412 1.00897216 1.02255316]]
```

The function `compute_dosage` also accepts the argument `ref` from which the reference
alleles can be specified. (Consult `help(bgen_reader.compute_dosage)` for the full
specification.)

Another example now querying specific locus and sample.

```python
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
```

```bash
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
```

## Troubleshooting

### fatal error: bgen.h: No such file or directory

This means that bgen C library is not installed (or could not be found). Please,
follow the instructions in https://github.com/limix/bgen to install it, and try
installing bgen-reader again.

## Problems

If you encounter any issue, please, [submit it](https://github.com/limix/bgen-reader-py/issues/new).

## Authors

* [Danilo Horta](https://github.com/horta)

## License

This project is licensed under the [MIT License](https://raw.githubusercontent.com/limix/bgen-reader-py/master/LICENSE.md).
