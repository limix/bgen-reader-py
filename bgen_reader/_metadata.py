from os.path import exists
from pathlib import Path

from ._bgen_file import bgen_file
from ._file import assert_file_exist, assert_file_readable
from ._string import make_sure_bytes


def create_metafile(bgen_filepath: Path, metafile_filepath: Path, verbose=True):
    r"""Create variants metadata file.

    Variants metadata file helps speed up subsequent reads of the associated
    bgen file.

    Parameters
    ----------
    bgen_filepath : str
        Bgen file path.
    metafile_file : str
        Metafile file path.
    verbose : bool
        ``True`` to show progress; ``False`` otherwise.

    Examples
    --------
    .. doctest::

        >>> import os
        >>> from bgen_reader import create_metafile, example_files
        >>>
        >>> with example_files("example.32bits.bgen") as filepath:
        ...     folder = os.path.dirname(filepath)
        ...     metafile_filepath = os.path.join(folder, filepath + ".metadata")
        ...
        ...     try:
        ...         create_metafile(filepath, metafile_filepath, verbose=False)
        ...     finally:
        ...         if os.path.exists(metafile_filepath):
        ...             os.remove(metafile_filepath)
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
        bgen.create_metafile(metafile_filepath, verbose)
