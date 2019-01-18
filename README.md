# bgen-reader

[![Travis](https://img.shields.io/travis/limix/bgen-reader-py.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.org/limix/bgen-reader-py) [![AppVeyor](https://img.shields.io/appveyor/ci/Horta/bgen-reader-py.svg?style=flat-square&label=windows%20build)](https://ci.appveyor.com/project/Horta/bgen-reader-py) [![Doc](https://readthedocs.org/projects/bgen-reader/badge/?version=latest&style=flat-square)](https://bgen-reader.readthedocs.io)

A [bgen file format](http://www.well.ox.ac.uk/~gav/bgen_format/) reader.

Bgen is a file format for storing large genetic datasets.
It supports both unphased genotypes and phased haplotype data with variable
ploidy and number of alleles. It was designed to provides a compact data
representation without sacrificing variant access performance.

This python package is a wrapper around the [bgen library](https://github.com/limix/bgen),
a low-memory footprint reader that efficiently reads bgen files.
It fully supports the bgen format specifications: 1.2 and 1.3;
as well as their optional compressed formats.

## Install

It can be installed vi [conda](https://conda.io/docs/)

```bash
conda install -c conda-forge bgen-reader
```

Or via pip

```bash
pip install bgen-reader
```

## Documentation

Please, refer to [https://bgen-reader.readthedocs.io](https://bgen-reader.readthedocs.io).


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
