typedef unsigned char byte;
typedef int_fast64_t  inti;
typedef double        real;

typedef struct string {
    inti  len;
    byte *str;
} string;

typedef struct Variant
{
    string  id;
    string  rsid;
    string  chrom;
    inti    position;
    inti    nalleles;
    string *allele_ids;
} Variant;

typedef struct VariantGenotype VariantGenotype;
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

VariantGenotype* open_variant_genotype(VariantIndexing *index,
                                       inti             variant_idx);

void             read_variant_genotype(VariantIndexing *index,
                                       VariantGenotype *vg,
                                       real            *probabilities);

inti   get_nalleles(VariantGenotype *vg);
inti   get_ploidy(VariantGenotype *vg);
inti   get_ncombs(VariantGenotype *vg);

void   close_variant_genotype(VariantIndexing *index,
                              VariantGenotype *vg);

void   free(void *);

string string_duplicate(const string s);

inti   sample_ids_presence(BGenFile *bgen);
