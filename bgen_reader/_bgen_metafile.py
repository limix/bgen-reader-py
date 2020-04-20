from collections import OrderedDict
from pathlib import Path
from threading import RLock

import dask.dataframe as dd
from cachetools import LRUCache, cached
from dask.delayed import delayed
from numpy import empty, uint16, uint32, uint64, zeros
from pandas import DataFrame

from ._ffi import ffi, lib


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

    @property
    def partition_size(self) -> int:
        return _ceildiv(self.nvariants, self.npartitions)

    def _inner_read_partition(self, index: int):
        # start = time()
        partition = lib.bgen_metafile_read_partition(self._bgen_metafile, index)
        # print(f"Elapsed: {time() - start} for bgen_metafile_read_partition")
        if partition == ffi.NULL:
            raise RuntimeError(f"Could not read partition {partition}.")

        nvariants = lib.bgen_partition_nvariants(partition)

        # start = time()
        position = empty(nvariants, dtype=uint32)
        nalleles = empty(nvariants, dtype=uint16)
        offset = empty(nvariants, dtype=uint64)
        vid_max_len = ffi.new("uint32_t[]", 1)
        rsid_max_len = ffi.new("uint32_t[]", 1)
        chrom_max_len = ffi.new("uint32_t[]", 1)
        allele_ids_max_len = ffi.new("uint32_t[]", 1)
        # print(f"Elapsed: {time() - start} empty")

        # start = time()
        position_ptr = ffi.cast("uint32_t *", ffi.from_buffer(position))
        nalleles_ptr = ffi.cast("uint16_t *", ffi.from_buffer(nalleles))
        offset_ptr = ffi.cast("uint64_t *", ffi.from_buffer(offset))
        lib.read_partition_part1(
            partition,
            position_ptr,
            nalleles_ptr,
            offset_ptr,
            vid_max_len,
            rsid_max_len,
            chrom_max_len,
            allele_ids_max_len,
        )
        # print(f"Elapsed: {time() - start} read_partition")

        # start = time()
        vid = zeros(nvariants, dtype=f"S{vid_max_len[0]}")
        rsid = zeros(nvariants, dtype=f"S{rsid_max_len[0]}")
        chrom = zeros(nvariants, dtype=f"S{chrom_max_len[0]}")
        allele_ids = zeros(nvariants, dtype=f"S{allele_ids_max_len[0]}")
        # print(f"Elapsed: {time() - start} create_strings")

        # start = time()
        lib.read_partition_part2(
            partition,
            ffi.from_buffer("char[]", vid),
            vid_max_len[0],
            ffi.from_buffer("char[]", rsid),
            rsid_max_len[0],
            ffi.from_buffer("char[]", chrom),
            chrom_max_len[0],
            ffi.from_buffer("char[]", allele_ids),
            allele_ids_max_len[0],
        )
        # print(f"Elapsed: {time() - start} read_partition2")
        lib.bgen_partition_destroy(partition)

        return nvariants, vid, rsid, chrom, position, nalleles, allele_ids, offset

    def read_partition(self, index: int):
        (
            nvariants,
            vid,
            rsid,
            chrom,
            position,
            nalleles,
            allele_ids,
            offset,
        ) = self._inner_read_partition(index)
        # start = time()
        data = OrderedDict(
            [
                ("id", vid.astype(str)),
                ("rsid", rsid.astype(str)),
                ("chrom", chrom.astype(str)),
                ("pos", position),
                ("nalleles", nalleles),
                ("allele_ids", allele_ids.astype(str)),
                ("vaddr", offset),
            ]
        )
        # print(f"Elapsed: {time() - start} for building OrderedDict")

        # start = time()
        df = DataFrame(data)
        # print(f"Elapsed: {time() - start} for building dataframe")

        # start = time()
        index_offset = self.partition_size * index
        df.index = range(index_offset, index_offset + nvariants)
        # print(f"Elapsed: {time() - start} for final arrangements")
        return df

    def create_variants(self):
        nvariants = self.nvariants
        npartitions = self.npartitions
        dfs = []
        index_base = 0
        part_size = self.partition_size
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


def _ceildiv(a, b):
    return -(-a // b)
