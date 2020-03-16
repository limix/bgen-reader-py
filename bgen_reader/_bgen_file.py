from contextlib import contextmanager
from math import floor, sqrt
from pathlib import Path

from pandas import Series
from numpy import asarray, full, nan, float64

from ._ffi import ffi, lib
from ._string import create_string, make_sure_bytes


class bgen_file:
    def __init__(self, filepath: Path):
        self._filepath = filepath
        self._bgen_file = None

    @property
    def filepath(self) -> Path:
        return self._filepath

    @property
    def nvariants(self) -> int:
        return lib.bgen_file_nvariants(self._bgen_file)

    @property
    def nsamples(self) -> int:
        return lib.bgen_file_nsamples(self._bgen_file)

    @property
    def contain_samples(self) -> bool:
        return lib.bgen_file_contain_samples(self._bgen_file)

    def read_samples(self, verbose: bool) -> Series:
        nsamples = self.nsamples
        bgen_samples = lib.bgen_file_read_samples(self._bgen_file, verbose)
        if bgen_samples == ffi.NULL:
            raise RuntimeError(f"Could not fetch samples from the bgen file.")

        try:
            ids = [
                create_string(lib.bgen_samples_get(bgen_samples, i))
                for i in range(nsamples)
            ]
        finally:
            lib.bgen_samples_destroy(bgen_samples)

        return Series(ids, dtype=str, name="id")

    def create_metafile(self, filepath: Path, verbose: bool):
        n = _estimate_best_npartitions(self.nvariants)

        metafile = lib.bgen_metafile_create(
            self._bgen_file, bytes(filepath), n, verbose
        )
        if metafile == ffi.NULL:
            raise RuntimeError(f"Error while creating metafile: {filepath}.")

        lib.bgen_metafile_close(metafile)

    def read_genotype(self, offset: int):
        genotype = lib.bgen_file_open_genotype(self._bgen_file, offset)
        if genotype == ffi.NULL:
            raise RuntimeError(f"Could not open genotype (offset {offset})")

        nsamples = self.nsamples
        ncombs = lib.bgen_genotype_ncombs(genotype)
        probs = full((nsamples, ncombs), nan, dtype=float64)
        lib.bgen_genotype_read(genotype, ffi.cast("double *", probs.ctypes.data))

        phased = lib.bgen_genotype_phased(genotype)
        ploidy = asarray(
            [lib.bgen_genotype_ploidy(genotype, i) for i in range(nsamples)], int
        )
        missing = asarray(
            [lib.bgen_genotype_missing(genotype, i) for i in range(nsamples)], bool
        )
        lib.bgen_genotype_close(genotype)
        return {"probs": probs, "phased": phased, "ploidy": ploidy, "missing": missing}

    def close(self):
        self.__exit__()

    def __enter__(self):
        self._bgen_file = lib.bgen_file_open(bytes(self._filepath))
        if self._bgen_file == ffi.NULL:
            raise RuntimeError(f"Could not open {self._filepath}.")
        return self

    def __exit__(self, *_):
        if self._bgen_file is not None:
            lib.bgen_file_close(self._bgen_file)


def _estimate_best_npartitions(nvariants: int):
    min_variants = 128
    m = max(min(min_variants, nvariants), floor(sqrt(nvariants)))
    return nvariants // m
