typedef struct bgen_string {
  int len;
  char *str;
} bgen_string;

struct bgen_var {
  bgen_string id;
  bgen_string rsid;
  bgen_string chrom;
  int position;
  int nalleles;
  bgen_string *allele_ids;
};

struct bgen_file *open_bgen(const char *filepath);

void close_bgen(struct bgen_file *bgen);

int get_nsamples(struct bgen_file *bgen);

int get_nvariants(struct bgen_file *bgen);

bgen_string *read_samples(struct bgen_file *bgen);

void free_samples(const struct bgen_file *bgen, bgen_string *samples);

struct bgen_var *read_variants(struct bgen_file *bgen, struct bgen_vi **index);

void free_variants(const struct bgen_file *bgen, struct bgen_var *variants);

void free_index(struct bgen_vi *index);

struct bgen_vg *open_variant_genotype(struct bgen_vi *index,
                                      size_t variant_idx);

void read_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg,
                           double *probabilities);

int get_nalleles(struct bgen_vg *vg);
int get_ploidy(struct bgen_vg *vg);
int get_ncombs(struct bgen_vg *vg);

void close_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg);

void free(void *);

bgen_string string_duplicate(const bgen_string s);

int sample_ids_presence(struct bgen_file *bgen);
