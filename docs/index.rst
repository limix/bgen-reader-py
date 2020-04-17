###########################
Bgen-reader's documentation
###########################

:Date: |today|
:Version: |version|

|Bgen| is a file format for storing large genetic datasets.
It supports both unphased genotypes and phased haplotype data with variable
ploidy and number of alleles. It was designed to provide a compact data
representation without sacrificing variant access performance.

This Python package is a wrapper around the |bgen library|,
a low-memory footprint reader that efficiently reads bgen files.
It fully supports the bgen format specifications: 1.2 and 1.3;
as well as their optional compressed formats.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   daskapi
   numpyapi
   env_variables

Comments and bugs
=================

You can get the source and open issues |on Github|.

.. |on Github| raw:: html

   <a href="https://github.com/limix/bgen-reader-py" target="_blank">on Github⧉</a>

.. |Bgen| raw:: html

   <a href="http://www.well.ox.ac.uk/~gav/bgen_format/" target="_blank">Bgen⧉</a>

.. |bgen library| raw:: html

   <a href="https://github.com/limix/bgen" target="_blank">bgen library⧉</a>
