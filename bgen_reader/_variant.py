from pathlib import Path
from threading import RLock

import dask.dataframe as dd
from cachetools import LRUCache, cached
from dask.delayed import delayed

from ._bgen_metafile import bgen_metafile


def create_variants(nvariants: int, metafile_filepath: Path):

    with bgen_metafile(metafile_filepath) as mf:
        npartitions = mf.npartitions

    dfs = []
    index_base = 0
    part_size = _get_partition_size(nvariants, npartitions)
    divisions = []
    for i in range(npartitions):
        divisions.append(index_base)
        d = delayed(read_partition)(metafile_filepath, i, index_base)
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
def read_partition(metafile_filepath, part, index_base):
    with bgen_metafile(metafile_filepath) as mf:
        return mf.read_partition(part, index_base)


def _get_partition_size(nvariants: int, npartitions: int):
    return _ceildiv(nvariants, npartitions)


def _ceildiv(a, b):
    return -(-a // b)
