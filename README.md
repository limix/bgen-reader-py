# Bgen-reader

[![Travis](https://img.shields.io/travis/limix/bgen-reader-py.svg?style=flat-square)](https://travis-ci.org/limix/bgen-reader-py)

A [BGEN file format](http://www.well.ox.ac.uk/~gav/bgen_format/) reader.

BGEN is a file format for storing large genetic datasets. It supports both unphased genotypes and phased haplotype data with variable ploidy and number of alleles. It was designed to provides a compact data representation without sacrificing variant access performance.

This Python package is a wrapper around the [bgen library](https://github.com/limix/bgen),
a low-memory footprint reader that efficiently reads BGEN files.
It fully supports all the BGEN format specifications: 1.1, 1.2, and 1.3;
as well as their optional compressed formats.

## Install

The recommended way of installing it is via
[conda](http://conda.pydata.org/docs/index.html)

```bash
conda install -c conda-forge bgen-reader
```

An alternative way would be via pip

```bash
pip install bgen-reader
```

In this case, you need to make sure you have both
[Zstandard](https://github.com/facebook/zstd) and [bgen](https://github.com/limix/bgen)
libraries properly installed.

## Usage

It is as simple as

```python
# example.py file
from bgen_reader import read_bgen

bgen = read_bgen("example.bgen", verbose=False)

print(bgen['variants'].head())
print(bgen['samples'].head())
print(len(bgen['genotype']))
print(bgen['genotype'][0].compute())
```

The output should something similar to

```
chrom       id  nalleles   pos    rsid
0    01  SNPID_2         2  2000  RSID_2
1    01  SNPID_3         2  3000  RSID_3
2    01  SNPID_4         2  4000  RSID_4
3    01  SNPID_5         2  5000  RSID_5
4    01  SNPID_6         2  6000  RSID_6
         id
0  sample_001
1  sample_002
2  sample_003
3  sample_004
4  sample_005
199
[[        nan         nan         nan]
[ 0.02780236  0.00863674  0.9635609 ]
[ 0.01736504  0.04968414  0.93295083]
...,
[ 0.01419069  0.02810669  0.95770262]
[ 0.91949463  0.05206298  0.02844239]
[ 0.00244141  0.98410029  0.0134583 ]]
```

## Problems

If you encounter any issue, please, [submit it](https://github.com/limix/bgen-reader-py/issues).

## Authors

* **Danilo Horta** - [https://github.com/horta](https://github.com/horta)

## License

This project is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details
