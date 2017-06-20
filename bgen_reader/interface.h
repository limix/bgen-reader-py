typedef struct PyArrayObject PyArrayObject;

typedef unsigned char byte;
typedef int_fast64_t  inti;
typedef double        real;

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
    inti  ncombs;
    real *probabilities;
} VariantGenotype;

typedef struct VariantIndexing VariantIndexing;
typedef struct BGenFile        BGenFile;


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

void           free_variant_genotypes(VariantGenotype *vg,
                                      inti             nvariants);

PyArrayObject* read_variant_genotype(VariantIndexing *indexing,
                                     inti             nsamples,
                                     inti             variant_idx);

void   free(void *);

string string_duplicate(const string s);

inti   sample_ids_presence(BGenFile *bgen);
