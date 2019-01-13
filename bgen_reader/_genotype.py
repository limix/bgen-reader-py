from threading import RLock

import dask
import dask.bag
from cachetools import LRUCache, cached
from numpy import asarray, float64, full, nan
from tqdm import trange

from ._bgen import bgen_file, bgen_metafile
from ._ffi import ffi, lib
from ._partition import read_partition


def map_genotype(bgen_filepath, metafile_filepath, verbose):
    with bgen_file(bgen_filepath) as bgen:
        nvariants = lib.bgen_nvariants(bgen)

    rg = _get_read_genotype(bgen_filepath, metafile_filepath)

    desc = "Mapping variants"
    return [
        rg(i, dask_key_name=str(i))
        for i in trange(nvariants, desc=desc, disable=not verbose)
    ]


def _get_read_genotype(bgen_filepath, metafile_filepath):
    @dask.delayed(nout=0, traverse=False, name="_read_genotype")
    def _read_genotype(i):

        with bgen_file(bgen_filepath) as bgen:
            nsamples = lib.bgen_nsamples(bgen)
            nvariants = lib.bgen_nvariants(bgen)
        with bgen_metafile(metafile_filepath) as mf:
            nparts = lib.bgen_metafile_nparts(mf)

        part_size = nvariants // nparts
        part = i // part_size
        index_base = (i % part_size) * part_size
        j = i % part_size
        p = read_partition(bgen_filepath, metafile_filepath, part, index_base)
        nsub_parts = _estimate_best_nsub_parts(nsamples, part_size)
        spart_size = max(1, part_size // nsub_parts)
        sub_part = j // spart_size
        m = j % spart_size
        start = sub_part * spart_size
        end = min(len(p), (sub_part + 1) * spart_size)
        vaddrs = tuple(p.iloc[start:end]["vaddr"].tolist())
        g = read_genotype_partition(
            bgen_filepath, metafile_filepath, vaddrs, sub_part, spart_size
        )
        return g[m]

    return _read_genotype


cache = LRUCache(maxsize=3)
lock = RLock()


@cached(cache, lock=lock)
def read_genotype_partition(
    bgen_filepath, metafile_filepath, vaddrs, sub_part, spart_size
):
    genotypes = []
    for vaddr in vaddrs:
        with bgen_file(bgen_filepath) as bgen:
            nsamples = lib.bgen_nsamples(bgen)
            vg = lib.bgen_open_genotype(bgen, vaddr)
            ncombs = lib.bgen_ncombs(vg)
            p = full((nsamples, ncombs), nan, dtype=float64)
            lib.bgen_read_genotype(bgen, vg, ffi.cast("double *", p.ctypes.data))
            phased = lib.bgen_phased(vg)
            ploidy = asarray([lib.bgen_ploidy(vg, i) for i in range(nsamples)], int)
            missing = asarray([lib.bgen_missing(vg, i) for i in range(nsamples)], bool)
            lib.bgen_close_genotype(vg)
            genotypes.append(
                {"probs": p, "phased": phased, "ploidy": ploidy, "missing": missing}
            )
    return genotypes


def _estimate_best_nsub_parts(nsamples, part_size):
    # Assume ideal block size, `bs`: 256KB
    # Assume 16 bytes per genotype per sample, `vs`
    # ideal nvariants to read: iv = bs / (vs * nsamples)
    # We then use iv to figure out in how many parts a partition will be subdivided
    # Let part_size be the number of variants in a partition
    # nsub_parts = min(int(part_size / iv), 1)
    bs = 256 * 1024
    vs = 16
    iv = bs / (vs * nsamples)
    return max(int(part_size / iv), 1)


# "probs": array([]),
# "phased": 0,
# "ploidy": array([2, 1, 2]),
# "missing": array([0, 0, 0]),
# }

# cache = LRUCache(maxsize=2)
# lock = RLock()

# genotype[5] = {
#   "probs": array([]),
#   "phased": 0,
#   "ploidy": array([2, 1, 2]),
#   "missing": array([0, 0, 0]),
# }

# ncombs = bgen_ncombs(vg)
#         max_ncombs = max(max_ncombs, ncombs)
#         ncombss.append(ncombs)

#         phased.append([bgen_phased(vg)] * nsamples)

#         ploidy.append([bgen_ploidy(vg, j) for j in range(nsamples)])

# missing.append([bgen_missing(vg, j) for j in range(nsamples)])

# @cached(cache, lock=lock)
# def _genotype_block(indexing, nsamples, variant_idx, nvariants):

#     max_ncombs = -inf
#     ncombss = []
#     variants = []
#     phased = []
#     ploidy = []
#     missing = []

#     for i in range(variant_idx, variant_idx + nvariants):
#         vg = bgen_open_variant_genotype(indexing[0], i)

#         ncombs = bgen_ncombs(vg)
#         max_ncombs = max(max_ncombs, ncombs)
#         ncombss.append(ncombs)

#         phased.append([bgen_phased(vg)] * nsamples)

#         ploidy.append([bgen_ploidy(vg, j) for j in range(nsamples)])

#         missing.append([bgen_missing(vg, j) for j in range(nsamples)])

#         g = full((nsamples, ncombs), nan, dtype=float64)

#         pg = ffi.cast("double *", g.ctypes.data)
#         bgen_read_variant_genotype(indexing[0], vg, pg)

#         bgen_close_variant_genotype(indexing[0], vg)

#         variants.append(g)

#     G = full((nvariants, nsamples, max_ncombs), nan, dtype=float64)

#     for i in range(0, nvariants):
#         G[i, :, : ncombss[i]] = variants[i]

#     phased = asarray(phased, int)
#     ploidy = asarray(ploidy, int)
#     missing = asarray(missing, int)

#     variant_idxs = range(variant_idx, variant_idx + nvariants)

#     data = stack([phased, ploidy, missing], axis=2)

#     coords = {
#         "variant": variant_idxs,
#         "sample": range(nsamples),
#         "data": ["phased", "ploidy", "missing"],
#     }
#     dims = ("variant", "sample", "data")
#     X = xr.DataArray(data, coords=coords, dims=dims)

#     return G, X
