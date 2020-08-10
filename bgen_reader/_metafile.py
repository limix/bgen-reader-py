import warnings
from pathlib import Path
from typing import Union

from ._bgen_file import bgen_file
from ._environment import BGEN_READER_CACHE_HOME
from ._file import (
    assert_file_exist,
    assert_file_readable,
    is_file_writable,
    path_to_filename,
)


def create_metafile(
    bgen_filepath: Union[str, Path],
    metafile_filepath: Union[str, Path],
    verbose: bool = True,
):
    """
    Create variants metadata file.

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


_metafile_nowrite_dir = """\
You don't have permission to write `{filepath}`. This might prevent speeding-up the reading process
in future runs.
"""


def infer_metafile_filepath(bgen_filepath: Path, suffix: str = ".metafile") -> Path:
    """
    Infer metafile filepath.

    The resulting file name will the file name of ``bgen_filepath`` with the appended ``suffix``.
    The root directory of the resulting filepath will be the directory of ``bgen_filepath`` if
    the user has appropriate permissions. It falls back to the directory

        BGEN_READER_CACHE_HOME / "metafile"

    if necessary.
    """
    metafile = bgen_filepath.with_suffix(bgen_filepath.suffix + suffix)
    if metafile.exists():
        try:
            assert_file_readable(metafile)
            return metafile
        except RuntimeError as e:
            warnings.warn(str(e), UserWarning)
            return BGEN_READER_CACHE_HOME / "metafile" / path_to_filename(metafile)
    else:
        if is_file_writable(metafile):
            return metafile

        warnings.warn(_metafile_nowrite_dir.format(filepath=metafile), UserWarning)
        return BGEN_READER_CACHE_HOME / "metafile" / path_to_filename(metafile)
