# bgen-reader

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

It can be installed via [conda](https://conda.io/docs/)

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

This project is licensed under the [MIT License](https://raw.githubusercontent.com/limix/bgen-reader-py/main/LICENSE.md).
