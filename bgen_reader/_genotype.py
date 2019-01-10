# "probs": array([]),
# "phased": 0,
# "ploidy": array([2, 1, 2]),
# "missing": array([0, 0, 0]),
# }

# cache = LRUCache(maxsize=2)
# lock = RLock()


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
