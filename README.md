
Bgen-reader
===========

|PyPI-Status| |Build-Status| |Codacy-Grade| |License-Badge|

A `BGEN file format`_ reader.

BGEN is a file format for storing large genetic datasets.
It supports both unphased genotypes and phased haplotype data with variable
ploidy and number of alleles. It was designed to provides a compact data
representation without sacrificing variant access performance.

This Python package is a wrapper around the `bgen library`_,
a low-memory footprint reader that efficiently reads BGEN files.
It fully supports all the BGEN format specifications: 1.1, 1.2, and 1.3;
as well as their optional compressed formats.

Install
-------

The recommended way of installing it is via conda_

.. code:: bash

    conda install -c conda-forge bgen-reader

An alternative way would be via pip_

.. code:: bash

    pip install bgen-reader

In this case, you need to make sure you have both Zstandard_ and bgen_
libraries properly installed.

Usage
-----

It is as simple as

.. code:: python

    # example.py file
    from bgen_reader import read_bgen

    bgen = read_bgen("example.bgen", verbose=False)

    print(bgen['variants'].head())
    print(bgen['samples'].head())
    print(len(bgen['genotype']))
    print(bgen['genotype'][0].compute())

The output should something similar to

.. code::

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

Problems
--------

If you encounter any issue, please, `submit it`_.

Authors
-------

* `Danilo Horta`_

License
-------

This project is licensed under the MIT License - see the `license file`_ for
details.

.. |Build-Status| image:: https://travis-ci.org/limix/bgen-reader-py.svg?branch=master
    :target: https://travis-ci.org/limix/bgen-reader-py

.. |Codacy-Grade| image:: https://api.codacy.com/project/badge/Grade/afb406c08b704f8a8722d8fe8e1b66f4
    :target: https://www.codacy.com/app/danilo.horta/bgen-reader-py?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=limix/bgen-reader-py&amp;utm_campaign=Badge_Grade

.. |PyPI-Status| image:: https://img.shields.io/pypi/v/bgen-reader.svg
    :target: https://pypi.python.org/pypi/bgen-reader

.. |License-Badge| image:: https://img.shields.io/pypi/l/bgen-reader.svg
    :target: https://raw.githubusercontent.com/bgen-reader-py/bgen-reader-py/master/LICENSE.txt

.. _License file: https://raw.githubusercontent.com/limix/bgen-reader-py/master/LICENSE.txt

.. _Danilo Horta: https://github.com/horta

.. _conda: http://conda.pydata.org/docs/index.html

.. _pip: https://pypi.python.org/pypi/pip

.. _Zstandard: https://github.com/facebook/zstd

.. _bgen: https://github.com/limix/bgen

.. _submit it: https://github.com/limix/bgen-reader-py/issues

.. _BGEN file format: http://www.well.ox.ac.uk/~gav/bgen_format/

.. _bgen library: https://github.com/limix/bgen
