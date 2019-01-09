from math import floor, sqrt
from os import W_OK, access
from os.path import abspath, dirname, exists

from ._ffi import ffi
from ._ffi.lib import (
    bgen_close,
    bgen_close_metafile,
    bgen_create_metafile,
    bgen_open,
    bgen_open_metafile,
    bgen_nvariants
)
from ._file import assert_file_exist, assert_file_readable
from ._misc import make_sure_bytes


def try_read_variants_metadata_file(bfile, mfilepath, index, v):
    if exists(mfilepath):
        variants = bgen_load_variants_metadata(bfile, mfilepath, index, v)
        if variants == ffi.NULL:
            if v == 1:
                msg = "Warning: could not read variants"
                msg += " metadata from {}.".format(mfilepath)
                print(msg)
            variants = bgen_read_variants_metadata(bfile, index, v)
    else:
        variants = bgen_read_variants_metadata(bfile, index, v)

    if variants == ffi.NULL:
        raise RuntimeError("Could not read variants metadata.")

    if not exists(mfilepath):
        if access(abspath(dirname(mfilepath)), W_OK):
            e = bgen_store_variants_metadata(bfile, variants, index[0], mfilepath)
            if e != 0 and v == 1:
                errmsg = "Warning: could not create"
                errmsg += " the metadata file {}.".format(abspath(mfilepath))
                print(errmsg)
        elif v == 1:
            errmsg = "Warning: you don't have permission to write"
            errmsg += " the metadata file {}.".format(abspath(mfilepath))
            print(errmsg)
    return variants


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

    bgen = bgen_open(make_sure_bytes(bgen_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not read {bgen_filepath}.")

    nparts = _estimate_best_npartitions(bgen_nvariants(bgen))
    metafile = bgen_create_metafile(bgen, metafile_filepath, nparts, verbose)
    if metafile == ffi.NULL:
        raise RuntimeError(f"Error while creating metafile: {metafile_filepath}.")

    if bgen_close_metafile(metafile) != 0:
        raise RuntimeError(f"Error while closing metafile: {metafile_filepath}.")

    bgen_close(bgen)


def _estimate_best_npartitions(nvariants):
    min_variants = 128
    m = max(min(min_variants, nvariants), floor(sqrt(nvariants)))
    return nvariants // m
