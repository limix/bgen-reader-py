from pathlib import Path
from typing import Union

from ._bgen_file import bgen_file
from ._file import assert_file_exist, assert_file_readable


def create_metafile(
    bgen_filepath: Union[str, Path],
    metafile_filepath: Union[str, Path],
    verbose: bool = True,
):
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
        >>> from bgen_reader import create_metafile, example_filepath
        >>>
        >>> filepath = example_filepath("example.32bits.bgen")
        >>> metafile_filepath = filepath.with_suffix(".metafile")
        >>>
        >>> try:
        ...     create_metafile(filepath, metafile_filepath, verbose=False)
        ... finally:
        ...     if metafile_filepath.exists():
        ...         os.remove(metafile_filepath)
    """
    bgen_filepath = Path(bgen_filepath)
    metafile_filepath = Path(metafile_filepath)

    assert_file_exist(bgen_filepath)
    assert_file_readable(bgen_filepath)

    if metafile_filepath.exists():
        raise ValueError(f"File {metafile_filepath} already exists.")

    with bgen_file(bgen_filepath) as bgen:
        bgen.create_metafile(metafile_filepath, verbose)
