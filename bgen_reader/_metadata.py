from math import floor, sqrt
from os.path import exists

from ._bgen import bgen_file
from ._ffi import ffi, lib
from ._file import assert_file_exist, assert_file_readable
from ._string import make_sure_bytes


def create_metafile(bgen_filepath, metafile_filepath, verbose=True):
    r"""Create variants metadata file.

    Variants metadata file helps speed up subsequent reads of the associated
    BGEN file.

    Parameters
    ----------
    bgen_filepath : str
        BGEN file path.
    metadata_file : str
        Metadata file path.
    verbose : bool
        ``True`` to show progress; ``False`` otherwise.
    """
    if verbose:
        verbose = 1
    else:
        verbose = 0

    bgen_filepath = make_sure_bytes(bgen_filepath)
    metafile_filepath = make_sure_bytes(metafile_filepath)

    assert_file_exist(bgen_filepath)
    assert_file_readable(bgen_filepath)

    if exists(metafile_filepath):
        raise ValueError(f"The file {metafile_filepath} already exists.")

    with bgen_file(bgen_filepath) as bgen:
        nparts = _estimate_best_npartitions(lib.bgen_nvariants(bgen))
        metafile = lib.bgen_create_metafile(bgen, metafile_filepath, nparts, verbose)
        if metafile == ffi.NULL:
            raise RuntimeError(f"Error while creating metafile: {metafile_filepath}.")

        if lib.bgen_close_metafile(metafile) != 0:
            raise RuntimeError(f"Error while closing metafile: {metafile_filepath}.")


def _estimate_best_npartitions(nvariants):
    min_variants = 128
    m = max(min(min_variants, nvariants), floor(sqrt(nvariants)))
    return nvariants // m
