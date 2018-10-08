from os import access, W_OK
from os.path import exists, dirname, abspath
from ._misc import make_sure_bytes, check_file_exist, check_file_readable

from ._ffi import ffi
from ._ffi.lib import (
    bgen_read_variants_metadata,
    bgen_create_variants_metadata_file,
    bgen_load_variants_metadata,
    bgen_store_variants_metadata,
)


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
        if access(dirname(mfilepath), W_OK):
            e = bgen_store_variants_metadata(
                bfile, variants, index[0], mfilepath
            )
            if e != 0 and v == 1:
                errmsg = "Warning: could not create"
                errmsg += " the metadata file {}.".format(abspath(mfilepath))
                print(errmsg)
        elif v == 1:
            errmsg = "Warning: you don't have permission to write"
            errmsg += " the metadata file {}.".format(abspath(mfilepath))
            print(errmsg)
    return variants


def create_metadata_file(bgen_filepath, metadata_filepath, verbose=True):
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
    metadata_filepath = make_sure_bytes(metadata_filepath)

    check_file_exist(bgen_filepath)
    check_file_readable(bgen_filepath)

    if exists(metadata_filepath):
        raise ValueError(
            "The file {} already exists.".format(metadata_filepath)
        )

    e = bgen_create_variants_metadata_file(
        bgen_filepath, metadata_filepath, verbose
    )

    if e != 0:
        raise RuntimeError("Error while creating metadata file: {}".format(e))
