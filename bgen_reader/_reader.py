import os
import warnings
from ._bgen import bgen_file
from ._ffi import lib
from ._file import (
    _get_temp_filepath,
    assert_file_exist,
    assert_file_exist2,
    assert_file_readable,
    is_file_writable,
    assert_file_readable2,
    BGEN_CACHE_HOME,
    path_to_filename,
    make_sure_dir_exist,
)
from ._genotype import create_genotypes
from ._metadata import create_metafile
from ._partition import create_variants
from ._samples import get_samples
from typing import Union
from pathlib import Path


def read_bgen(
    filepath: Union[str, Path],
    metafile_filepath=None,
    samples_filepath=None,
    verbose=True,
):
    r""" Read a given BGEN file.

    Parameters
    ----------
    filepath : str
        A bgen file path.
    metafile_filepath : str, optional
        If ``None``, it will try to read the ``filepath + ".metadata"`` file. If this is
        not possible, it will create one. It tries to create one at
        ``filepath + ".metadata"``. If that is also no possible, it tries to create one
        at a temporary folder.
    samples_filepath : str, optional
        A sample file in `gen format <https://goo.gl/bCzo7m>`_.
        If ``samples_filepath`` is provided, sample ids are read from this file.
        Otherwise, it reads from the bgen file itself if possible. Defaults to ``None``.
    verbose : bool, optional
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

        >>> from bgen_reader import example_files, read_bgen
        >>>
        >>> with example_files("haplotypes.bgen") as filepath:
        ...     bgen = read_bgen(filepath, verbose=False)
        ...     variants = bgen["variants"]
        ...     samples = bgen["samples"]
        ...
        ...     v = variants.loc[0].compute()
        ...     g = bgen["genotype"][0].compute()
        ...     print(v)
        ...     print(samples)
        ...     print(g["probs"][0])
             id rsid chrom  pos  nalleles allele_ids  vaddr
        0  SNP1  RS1     1    1         2        A,G    102
        0    sample_0
        1    sample_1
        2    sample_2
        3    sample_3
        Name: id, dtype: object
        [1. 0. 1. 0.]
    """

    filepath = Path(filepath)
    assert_file_exist2(filepath)
    assert_file_readable2(filepath)

    metafile_filepath = _get_valid_metafile_filepath(filepath, metafile_filepath)
    if not os.path.exists(metafile_filepath):
        if verbose:
            print(
                f"We will create the metafile `{metafile_filepath}`. This file will "
                "speed up further\nreads and only need to be created once. So, please, "
                "bear with me."
            )
        create_metafile(filepath, metafile_filepath, verbose)
    elif os.path.getmtime(metafile_filepath) < os.path.getmtime(filepath):
        if verbose:
            msg = f"File `{filepath}` has been modified after the creation of `{metafile_filepath}`.\n"
            msg += "We will therefore recreate the metadata file. So, please, bear with me."
            print(msg)
        os.unlink(metafile_filepath)
        create_metafile(filepath, metafile_filepath, verbose)

    with bgen_file(filepath) as bgen:
        nvariants = lib.bgen_file_nvariants(bgen)

    samples = get_samples(filepath, samples_filepath, verbose)
    variants = create_variants(nvariants, metafile_filepath)
    genotype = create_genotypes(filepath, metafile_filepath, verbose)

    return dict(variants=variants, samples=samples, genotype=genotype)


_metafile_not_found = """\
Metafile `{filepath}` does not exist. If you wish to create a metafile in a custom
location, please use `bgen_reader.create_metafile`.
"""

_metafile_nowrite_dir = """\
You don't have permission to write `{filepath}`. This might prevent speeding-up the reading process
in future runs.
"""


def _get_valid_metafile_filepath(bgen_filepath, metafile_filepath):
    if metafile_filepath is None:
        return _infer_metafile_filepath(bgen_filepath)
    else:
        try:
            assert_file_exist(metafile_filepath)
        except FileNotFoundError as e:
            fp = metafile_filepath
            warnings.warn(_metafile_not_found.format(filepath=fp), UserWarning)
            raise e
    return metafile_filepath


def _infer_metafile_filepath(bgen_filepath):
    metafile = bgen_filepath.with_suffix(".metadata")
    if metafile.exists():
        try:
            assert_file_readable2(metafile)
            return metafile
        except RuntimeError as e:
            warnings.warn(str(e), UserWarning)
            return BGEN_CACHE_HOME / path_to_filename(metafile)
    else:
        if is_file_writable(metafile):
            return metafile

        warnings.warn(_metafile_nowrite_dir.format(filepath=metafile), UserWarning)
        return BGEN_CACHE_HOME / path_to_filename(metafile)
