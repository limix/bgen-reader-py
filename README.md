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


## Troubleshooting

### fatal error: bgen.h: No such file or directory

This means that bgen C library is not installed (or could not be found). Please,
follow the instructions in <https://github.com/limix/bgen> to install it, and try
installing bgen-reader again.

## Problems

If you encounter any issue, please, [submit it](https://github.com/limix/bgen-reader-py/issues/new).

## Authors

* [Danilo Horta](https://github.com/horta)

## License

This project is licensed under the [MIT License](https://raw.githubusercontent.com/limix/bgen-reader-py/master/LICENSE.md).
