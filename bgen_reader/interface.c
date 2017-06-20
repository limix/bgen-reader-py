#include <stdlib.h>

#include "Python.h"
#include "bgen/bgen.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ndarrayobject.h"

BGenFile* open_bgen(const byte *filepath)
{
    return bgen_open(filepath);
}

void close_bgen(BGenFile *bgen)
{
    bgen_close(bgen);
}

inti get_nsamples(BGenFile *bgen)
{
    return bgen_nsamples(bgen);
}

inti get_nvariants(BGenFile *bgen)
{
    return bgen_nvariants(bgen);
}

string* read_samples(BGenFile *bgen)
{
    return bgen_read_samples(bgen);
}

void free_samples(const BGenFile *bgen,
                  string         *samples)
{
    bgen_free_samples(bgen, samples);
}

Variant* read_variants(BGenFile         *bgen,
                       VariantIndexing **index)
{
    return bgen_read_variants(bgen, index);
}

void free_variants(const BGenFile *bgen,
                   Variant        *variants)
{
    bgen_free_variants(bgen, variants);
}

void free_indexing(VariantIndexing *index)
{
    return bgen_free_indexing(index);
}

VariantGenotype* read_variant_genotypes(VariantIndexing *indexing,
                                        inti             variant_start,
                                        inti             variant_end)
{
    return bgen_read_variant_genotypes(indexing, variant_start, variant_end);
}

PyArrayObject* read_variant_genotype(VariantIndexing *indexing,
                                     inti             nsamples,
                                     inti             variant_idx)
{
    VariantGenotype *vg = bgen_read_variant_genotypes(indexing,
                                                      variant_idx,
                                                      variant_idx + 1);

    npy_intp dims[2]   = { nsamples, vg->ncombs };
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNewFromData(2,
                                                                    dims,
                                                                    NPY_FLOAT64,
                                                                    vg->probabilities);

    PyArray_ENABLEFLAGS(arr, NPY_ARRAY_OWNDATA);
    free(vg);

    return arr;
}

// ctypedef np.int32_t DTYPE_t
//
// cdef extern from "numpy/arrayobject.h":
//     void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
//
// cdef data_to_numpy_array_with_spec(void * ptr, np.npy_intp N, int t):
//     cdef np.ndarray[DTYPE_t, ndim=1] arr = np.PyArray_SimpleNewFromData(1,
// &N, t, ptr)
//     PyArray_ENABLEFLAGS(arr, np.NPY_OWNDATA)
//     return arr

void free_variant_genotypes(VariantGenotype *vg,
                            inti             nvariants)
{
    bgen_free_variant_genotypes(vg, nvariants);
}

string string_duplicate(string s)
{
    string r;

    r.str = malloc(s.len);
    memcpy(r.str, s.str, s.len);
    r.len = s.len;
    return r;
}

inti sample_ids_presence(BGenFile *bgen)
{
    return bgen_sample_ids_presence(bgen);
}
