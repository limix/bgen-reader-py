from math import floor, sqrt
from pathlib import Path

from numpy import float64, full, nan, uint8, zeros
from pandas import Series

from ._ffi import ffi, lib


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

    def read_samples(self) -> Series:
        nsamples = self.nsamples
        bgen_samples = lib.bgen_file_read_samples(self._bgen_file)
        if bgen_samples == ffi.NULL:
            raise RuntimeError("Could not fetch samples from the bgen file.")

        try:
            samples_max_len = ffi.new("uint32_t[]", 1)
            lib.read_samples_part1(bgen_samples, nsamples, samples_max_len)
            samples = zeros(nsamples, dtype=f"S{samples_max_len[0]}")
            lib.read_samples_part2(
                bgen_samples,
                nsamples,
                ffi.from_buffer("char[]", samples),
                samples_max_len[0],
            )
        finally:
            lib.bgen_samples_destroy(bgen_samples)

        return Series(samples, dtype=str, name="id")

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

        ploidy = full(nsamples, 0, dtype=uint8)
        lib.read_ploidy(genotype, ffi.cast("uint8_t *", ploidy.ctypes.data), nsamples)

        missing = full(nsamples, 0, dtype=bool)
        lib.read_missing(genotype, ffi.cast("bool *", missing.ctypes.data), nsamples)

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
