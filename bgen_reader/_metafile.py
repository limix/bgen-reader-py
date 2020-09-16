import warnings
from pathlib import Path
from threading import RLock
from typing import Union

import dask.dataframe as dd
from cachetools import LRUCache, cached
from cbgen import bgen_file, bgen_metafile
from cbgen.typing import Partition
from dask.delayed import delayed
from pandas import DataFrame

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


def create_variants(filepath: Path, nvariants: int, npartitions: int, part_size: int):
    dfs = []
    index_base = 0
    divisions = []
    for i in range(npartitions):
        divisions.append(index_base)
        d = delayed(_read_partition)(filepath, i)
        dfs.append(d)
        index_base += part_size
    divisions.append(nvariants - 1)
    meta = [
        ("id", str),
        ("rsid", str),
        ("chrom", str),
        ("pos", int),
        ("nalleles", int),
        ("allele_ids", str),
        ("vaddr", int),
    ]
    df = dd.from_delayed(dfs, meta=dd.utils.make_meta(meta), divisions=divisions)
    return df


cache = LRUCache(maxsize=3)
lock = RLock()


@cached(cache, lock=lock)
def _read_partition(filepath: Path, partition: int) -> DataFrame:
    with bgen_metafile(filepath) as mf:
        part: Partition = mf.read_partition(partition)
    v = part.variants
    data = {
        "id": v.id.astype(str),
        "rsid": v.rsid.astype(str),
        "chrom": v.chromosome.astype(str),
        "pos": v.position.astype(int),
        "nalleles": v.nalleles.astype(int),
        "allele_ids": v.allele_ids.astype(str),
        "vaddr": v.offset.astype(int),
    }
    df = DataFrame(data)
    return df[["id", "rsid", "chrom", "pos", "nalleles", "allele_ids", "vaddr"]]
