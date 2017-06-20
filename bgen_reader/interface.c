#include <stdlib.h>
#include <stdio.h>

#include "bgen/bgen.h"

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

// VariantGenotype* read_variant_genotypes(VariantIndexing *indexing,
//                                         inti             variant_start,
//                                         inti             variant_end)
// {
//     return bgen_read_variant_genotypes(indexing, variant_start, variant_end);
// }

// PyArrayObject* read_variant_genotype(VariantIndexing *indexing,
//                                      inti             nsamples,
//                                      inti             variant_idx)
// {
//     printf("Ponto 1\n"); fflush(stdout);
//     VariantGenotype *vg = bgen_read_variant_genotypes(indexing,
//                                                       variant_idx,
//                                                       variant_idx + 1);
//     printf("Ponto 2\n"); fflush(stdout);
//     npy_intp dims[2] = { nsamples, vg->ncombs };
//     printf("Ponto 3\n"); fflush(stdout);
//     printf("nsamples: %lld\n", nsamples); fflush(stdout);
//     printf("ncombs: %lld\n",   vg->ncombs); fflush(stdout);
//     PyObject *arr = PyArray_SimpleNewFromData(2,
//                                               dims,
//                                               NPY_FLOAT64,
//                                               vg->probabilities);
//
//     // PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNewFromData(2,
//     //                                                                 dims,
//     //
//     //
//     //
//                                                            NPY_FLOAT64,
//     //
//     //
//     //
//     //
//     //
//                                                      vg->probabilities);
//     printf("Ponto 4\n"); fflush(stdout);
//     PyArray_ENABLEFLAGS(arr, NPY_ARRAY_OWNDATA);
//     printf("Ponto 5\n"); fflush(stdout);
//     free(vg);
//     printf("Ponto 6\n"); fflush(stdout);
//
//     return arr;
// }

// void free_variant_genotypes(VariantGenotype *vg,
//                             inti             nvariants)
// {
//     bgen_free_variant_genotypes(vg, nvariants);
// }

VariantGenotype* open_variant_genotype(VariantIndexing *indexing,
                                       inti             variant_idx)
{
    return bgen_open_variant_genotype(indexing, variant_idx);
}

void read_variant_genotype(VariantIndexing *indexing,
                           VariantGenotype *vg,
                           real            *probabilities)
{
    bgen_read_variant_genotype(indexing, vg, probabilities);
}

inti variant_genotype_nsamples(VariantGenotype *vg)
{
    return bgen_variant_genotype_nsamples(vg);
}

inti variant_genotype_nalleles(VariantGenotype *vg)
{
    return bgen_variant_genotype_nalleles(vg);
}

inti variant_genotype_ploidy(VariantGenotype *vg)
{
    return bgen_variant_genotype_ploidy(vg);
}

inti variant_genotype_ncombs(VariantGenotype *vg)
{
    return bgen_variant_genotype_ncombs(vg);
}

void close_variant_genotype(VariantIndexing *indexing,
                            VariantGenotype *vg)
{
    bgen_close_variant_genotype(indexing, vg);
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
