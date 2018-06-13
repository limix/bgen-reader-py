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

### Unphased genotype

```python
>>> from bgen_reader import read_bgen
>>>
>>> bgen = read_bgen("example.bgen", verbose=False)
>>>
>>> print(bgen["variants"].head())
```
```
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
```
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
```
199
```
```python
>>> p = bgen["genotype"][0].compute()
>>> print(p)
```
```
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
```
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
```
     id rsid chrom  pos  nalleles allele_ids
0  SNP1  RS1     1    1         2        A,G
1  SNP2  RS2     1    2         2        A,G
2  SNP3  RS3     1    3         2        A,G
3  SNP4  RS4     1    4         2        A,G
```
```python
>>> print(bgen["samples"].head())
```
```
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
```
[0. 1. 1. 0.]
```
```python
>>> # Is it a phased one?
>>> print(bgen["X"][0, 1].compute().sel(data="phased").item())
```
```
1
```
```python
>>> # How many haplotypes?
>>> print(bgen["X"][0, 1].compute().sel(data="ploidy").item())
```
```
2
```
```python
>>> # And how many alleles?
>>> print(bgen["variants"].loc[0, "nalleles"])
```
```
2
```
```python
>>> # Therefore, the first haplotype has probability 100%
>>> # of having the allele
>>> print(bgen["variants"].loc[0, "allele_ids"].split(",")[1])
```
```
G
```
```python
>>> # And the second haplotype has probability 100% of having
>>> # the first allele
>>> print(bgen["variants"].loc[0, "allele_ids"].split(",")[0])
```
```
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
```
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
```
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
```
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
```
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
```
[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
```
```python
>>> # of size
>>> print(len(p))
```
```
36
```
```python
>>> # Since the 9th variant for the 4th individual is
>>> # unphased,
>>> print(bgen["X"][8, 3].compute().sel(data="phased").item())
```
```
0
```
```python
>>> # the estimated probabilities imply the dosage
>>> # (or expected number of alleles)
>>> print(convert_to_dosage(p, nalleles, ploidy))
```
```
[0. 1. 0. 0. 0. 1. 0. 0.]
```

## Problems

If you encounter any issue, please, [submit it](https://github.com/limix/bgen-reader-py/issues/new).

## Authors

* [Danilo Horta](https://github.com/horta)

## License

This project is licensed under the [MIT License](https://raw.githubusercontent.com/limix/bgen-reader-py/master/LICENSE.md).
