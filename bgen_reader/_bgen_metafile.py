from contextlib import contextmanager
from pathlib import Path
from pandas import DataFrame

from ._ffi import ffi, lib
from ._string import create_string


class bgen_metafile2:
    def __init__(self, filepath: Path):
        self._filepath = filepath
        self._bgen_metafile = None

    @property
    def npartitions(self) -> int:
        return lib.bgen_metafile_npartitions(self._bgen_metafile)

    @property
    def nvariants(self) -> int:
        return lib.bgen_metafile_nvariants(self._bgen_metafile)

    def read_partition(self, partition: int, index_base: int):
        partition = lib.bgen_metafile_read_partition(self._bgen_metafile, partition)
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
        df.index = range(index_base, index_base + nvariants)

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


def _read_allele_ids(allele_ids, nalleles):
    alleles = [create_string(allele_ids[i]) for i in range(nalleles)]
    return ",".join(alleles)
