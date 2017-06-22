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

VariantGenotype* open_variant_genotype(VariantIndexing *index,
                                       inti             variant_idx)
{
    return bgen_open_variant_genotype(index, variant_idx);
}

void read_variant_genotype(VariantIndexing *index,
                           VariantGenotype *vg,
                           real            *probabilities)
{
    bgen_read_variant_genotype(index, vg, probabilities);
}

inti get_nalleles(VariantGenotype *vg)
{
    return bgen_nalleles(vg);
}

inti get_ploidy(VariantGenotype *vg)
{
    return bgen_ploidy(vg);
}

inti get_ncombs(VariantGenotype *vg)
{
    return bgen_ncombs(vg);
}

void close_variant_genotype(VariantIndexing *index,
                            VariantGenotype *vg)
{
    bgen_close_variant_genotype(index, vg);
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
