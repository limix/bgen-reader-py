# Bgen-reader

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

```
pip install bgen-reader
```

In this case, you need to make sure you have both
[Zstandard](https://github.com/facebook/zstd) and [bgen](https://github.com/limix/bgen)
libraries properly installed.

## Problems

If you encounter any issue, please, [submit it](https://github.com/limix/bgen-reader-py/issues).

## Authors

* **Danilo Horta** - [https://github.com/horta](https://github.com/horta)

## License

This project is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details
