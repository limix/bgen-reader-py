typedef unsigned char byte;
typedef int_fast64_t  inti;
typedef double        real;

typedef struct string          string;
typedef struct BGenFile        BGenFile;
typedef struct Variant         Variant;
typedef struct VariantGenotype VariantGenotype;
typedef struct VariantIndexing VariantIndexing;

BGenFile       * open_bgen(const byte *filepath);

void             close_bgen(BGenFile *bgen);

inti             get_nsamples(BGenFile *bgen);

inti             get_nvariants(BGenFile *bgen);

string         * read_samples(BGenFile *bgen);

void             free_samples(const BGenFile *bgen,
                              string         *samples);

Variant        * read_variants(BGenFile         *bgen,
                               VariantIndexing **index);

void             free_variants(const BGenFile *bgen,
                               Variant        *variants);

void             free_indexing(VariantIndexing *index);

VariantGenotype* read_variant_genotypes(VariantIndexing *indexing,
                                        inti             variant_start,
                                        inti             variant_end);

void free_variant_genotypes(VariantGenotype *vg,
                            inti             nvariants);

void free(void *);
