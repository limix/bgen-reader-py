import os
from pathlib import Path
from typing import Optional, Union

from cbgen import bgen_file, bgen_metafile
from pandas import Series

from ._file import assert_file_exist, assert_file_readable
from ._genotype import create_genotypes
from ._metafile import create_metafile, create_variants, infer_metafile_filepath
from ._samples import generate_samples, read_samples_file


def read_bgen(
    filepath: Union[str, Path],
    metafile_filepath: Optional[Union[str, Path]] = None,
    samples_filepath: Optional[Union[str, Path]] = None,
    verbose: bool = True,
):
    """
    Read a given BGEN file.

    Parameters
    ----------
    filepath
        Bgen file path.
    metafile_filepath
        File path to the corresponding metafile. A metafile can be created by calling
        :func:`bgen_reader.create_metafile`. If ``None``, a metafile will be automatically created.
        Defaults to ``None``.
    samples_filepath
        Path to a `sample format`_ file or ``None`` to read samples from the bgen file itself.
        Defaults to ``None``.
    verbose
        ``True`` to show progress; ``False`` otherwise. Defaults to ``True``.

    Returns
    -------
    variants : :class:`dask.dataFrame.DataFrame`
        Variant position, chromosomes, rsids, etc.
    samples : :class:`pandas.Series`
        Sample identifications.
    genotype : list
        List of genotypes.

    Examples
    --------
    .. doctest::

        >>> from bgen_reader import example_filepath, read_bgen
        >>>
        >>> bgen = read_bgen(example_filepath("haplotypes.bgen"), verbose=False)
        >>> variants = bgen["variants"]
        >>> samples = bgen["samples"]
        >>>
        >>> v = variants.loc[0].compute()
        >>> g = bgen["genotype"][0].compute()
        >>> print(v)
             id rsid chrom  pos  nalleles allele_ids  vaddr
        0  SNP1  RS1     1    1         2        A,G    102
        >>> print(samples)
        0    sample_0
        1    sample_1
        2    sample_2
        3    sample_3
        Name: id, dtype: object
        >>> print(g["probs"][0])
        [1. 0. 1. 0.]

    .. _sample format: https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html
    """

    filepath = Path(filepath)
    assert_file_exist(filepath)
    assert_file_readable(filepath)

    if metafile_filepath is None:
        metafile_filepath = infer_metafile_filepath(filepath)
    else:
        metafile_filepath = Path(metafile_filepath)
        assert_file_exist(metafile_filepath)
        assert_file_readable(filepath)

    if not metafile_filepath.exists():
        if verbose:
            print(
                f"We will create the metafile `{metafile_filepath}`. This file will "
                "speed up further\nreads and only need to be created once. So, please, "
                "bear with me."
            )
        create_metafile(filepath, metafile_filepath, verbose)
    elif os.path.getmtime(metafile_filepath) < os.path.getmtime(filepath):
        from ._genotype import cache as bgencache
        from ._metafile import cache as metacache

        metacache.clear()
        bgencache.clear()

        if verbose:
            print(
                f"File `{filepath}` has been modified after the creation of `{metafile_filepath}`."
                "\nWe will therefore recreate the metadata file. So, please, bear with me."
            )
        os.unlink(metafile_filepath)
        create_metafile(filepath, metafile_filepath, verbose)

    with bgen_file(filepath) as bgen:
        samples = _get_samples(bgen, samples_filepath, verbose)

        with bgen_metafile(metafile_filepath) as mf:
            nvariants = mf.nvariants
            npartitions = mf.npartitions
            part_size = mf.partition_size
            variants = create_variants(
                metafile_filepath, nvariants, npartitions, part_size
            )

        genotype = create_genotypes(bgen, metafile_filepath, verbose)

    return dict(variants=variants, samples=samples, genotype=genotype)


def _get_samples(bgen, sample_file, verbose) -> Series:
    if sample_file is None:
        if bgen.contain_samples:
            samples = bgen.read_samples().astype(str)
            return Series(samples, dtype=str, name="id")
        else:
            return generate_samples(bgen.nsamples)
    else:
        samples_filepath = Path(sample_file)
        assert_file_exist(samples_filepath)
        assert_file_readable(samples_filepath)
        return read_samples_file(samples_filepath, verbose)
