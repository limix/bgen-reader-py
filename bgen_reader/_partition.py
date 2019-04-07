from threading import RLock

import dask.dataframe as dd
from cachetools import LRUCache, cached
from dask.delayed import delayed
from pandas import DataFrame

from ._bgen import bgen_file, bgen_metafile
from ._ffi import ffi, lib
from ._string import bgen_str_to_str


def map_metadata(bgen_filepath, metafile_filepath):
    with bgen_metafile(metafile_filepath) as mf:
        nparts = lib.bgen_metafile_npartitions(mf)
    with bgen_file(bgen_filepath) as bgen:
        nvariants = lib.bgen_nvariants(bgen)
    dfs = []
    index_base = 0
    part_size = get_partition_size(bgen_filepath, metafile_filepath)
    divisions = []
    for i in range(nparts):
        divisions.append(index_base)
        d = delayed(read_partition)(bgen_filepath, metafile_filepath, i, index_base)
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
def read_partition(bgen_filepath, metafile_filepath, part, index_base):
    with bgen_metafile(metafile_filepath) as mf:

        nvariants_ptr = ffi.new("int *")
        metadata = lib.bgen_read_partition(mf, part, nvariants_ptr)
        nvariants = nvariants_ptr[0]
        variants = []
        for i in range(nvariants):
            id_ = bgen_str_to_str(metadata[i].id)
            rsid = bgen_str_to_str(metadata[i].rsid)
            chrom = bgen_str_to_str(metadata[i].chrom)
            pos = metadata[i].position
            nalleles = metadata[i].nalleles
            allele_ids = _read_allele_ids(metadata[i])
            vaddr = metadata[i].vaddr
            variants.append([id_, rsid, chrom, pos, nalleles, allele_ids, vaddr])

        index = range(index_base, index_base + nvariants)
        variants = DataFrame(
            variants,
            index=index,
            columns=["id", "rsid", "chrom", "pos", "nalleles", "allele_ids", "vaddr"],
            dtype=str,
        )
        variants["pos"] = variants["pos"].astype(int)
        variants["nalleles"] = variants["nalleles"].astype(int)
        variants["vaddr"] = variants["vaddr"].astype(int)

    return variants


def get_partition_size(bgen_filepath, metafile_filepath):
    with bgen_file(bgen_filepath) as bgen:
        nvariants = lib.bgen_nvariants(bgen)
    with bgen_metafile(metafile_filepath) as mf:
        nparts = lib.bgen_metafile_npartitions(mf)
    return _ceildiv(nvariants, nparts)


def _ceildiv(a, b):
    return -(-a // b)


def _read_allele_ids(metadata):
    n = metadata.nalleles
    alleles = [bgen_str_to_str(metadata.allele_ids[i]) for i in range(n)]
    return ",".join(alleles)
