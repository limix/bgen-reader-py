from pathlib import Path

from pandas import DataFrame

from ._ffi import ffi, lib
from threading import RLock
from ._string import create_string
from dask.delayed import delayed
import dask.dataframe as dd
from cachetools import LRUCache, cached


class bgen_metafile:
    def __init__(self, filepath: Path):
        self._filepath = filepath
        self._bgen_metafile = None

    @property
    def filepath(self) -> Path:
        return self._filepath

    @property
    def npartitions(self) -> int:
        return lib.bgen_metafile_npartitions(self._bgen_metafile)

    @property
    def nvariants(self) -> int:
        return lib.bgen_metafile_nvariants(self._bgen_metafile)

    def read_partition(self, index: int):
        partition = lib.bgen_metafile_read_partition(self._bgen_metafile, index)
        if partition == ffi.NULL:
            raise RuntimeError(f"Could not read partition {partition}.")

        nvariants = self.nvariants
        variants = []
        for i in range(nvariants):
            variant = lib.bgen_partition_get_variant(partition, i)
            id_ = create_string(variant[0].id)
            rsid = create_string(variant[0].rsid)
            chrom = create_string(variant[0].chrom)
            pos = variant[0].position
            nalleles = variant[0].nalleles
            allele_ids = _read_allele_ids(variant[0].allele_ids, variant[0].nalleles)
            offset = variant[0].genotype_offset
            variants.append([id_, rsid, chrom, pos, nalleles, allele_ids, offset])

        df = DataFrame(
            variants,
            columns=["id", "rsid", "chrom", "pos", "nalleles", "allele_ids", "vaddr"],
            dtype=str,
        )
        df["pos"] = df["pos"].astype("uint32")
        df["nalleles"] = df["nalleles"].astype("uint16")
        df["vaddr"] = df["vaddr"].astype("uint64")

        part_size = _get_partition_size(nvariants, self.npartitions)
        index_offset = part_size * index
        df.index = range(index_offset, index_offset + nvariants)

        return df

    def create_variants(self):
        nvariants = self.nvariants
        npartitions = self.npartitions
        dfs = []
        index_base = 0
        part_size = _get_partition_size(nvariants, npartitions)
        divisions = []
        for i in range(npartitions):
            divisions.append(index_base)
            d = delayed(read_partition)(self._filepath, i)
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

    def close(self):
        self.__exit__()

    def __enter__(self):
        self._bgen_metafile = lib.bgen_metafile_open(bytes(self._filepath))
        if self._bgen_metafile == ffi.NULL:
            raise RuntimeError(f"Could not open {self._filepath}.")
        return self

    def __exit__(self, *_):
        if self._bgen_metafile is not None:
            lib.bgen_metafile_close(self._bgen_metafile)


cache = LRUCache(maxsize=3)
lock = RLock()


@cached(cache, lock=lock)
def read_partition(metafile_filepath: Path, partition: int):
    with bgen_metafile(metafile_filepath) as metafile:
        return metafile.read_partition(partition)


def _read_allele_ids(allele_ids, nalleles):
    alleles = [create_string(allele_ids[i]) for i in range(nalleles)]
    return ",".join(alleles)


def _get_partition_size(nvariants: int, npartitions: int):
    return _ceildiv(nvariants, npartitions)


def _ceildiv(a, b):
    return -(-a // b)
