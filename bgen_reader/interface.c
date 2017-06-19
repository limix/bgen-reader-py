#include <stdlib.h>

#include "bgen/bgen.h"

typedef struct string {
    inti  len;
    byte *str;
} string;

typedef struct BGenFile BGenFile;

typedef struct Variant
{
    string  id;
    string  rsid;
    string  chrom;
    inti    position;
    inti    nalleles;
    string *allele_ids;
} Variant;

typedef struct VariantGenotype
{
    inti  ploidy;
    real *probabilities;
} VariantGenotype;


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

void free_variant_genotypes(VariantGenotype *vg,
                            inti             nvariants)
{
    bgen_free_variant_genotypes(vg, nvariants);
}
