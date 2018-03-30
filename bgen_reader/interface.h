struct bgen_string {
  int len;
  char *str;
};

struct bgen_var {
  struct bgen_string id;
  struct bgen_string rsid;
  struct bgen_string chrom;
  int position;
  int nalleles;
  struct bgen_string *allele_ids;
};

struct bgen_file *open_bgen(const char *filepath);

void close_bgen(struct bgen_file *bgen);

int get_nsamples(const struct bgen_file *bgen);

int get_nvariants(const struct bgen_file *bgen);

int sample_ids_presence(const struct bgen_file *bgen);

struct bgen_string *read_samples(struct bgen_file *bgen, int verbose);

void free_samples(const struct bgen_file *bgen, struct bgen_string *samples);

struct bgen_var *read_variants_metadata(struct bgen_file *bgen,
                                        struct bgen_vi **index, int verbose);

void free_variants_metadata(const struct bgen_file *bgen,
                            struct bgen_var *variants);

void free_index(struct bgen_vi *index);

struct bgen_vg *open_variant_genotype(struct bgen_vi *index,
                                      size_t variant_idx);

void read_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg,
                           double *probabilities);

int get_nalleles(const struct bgen_vg *vg);
int get_ploidy(const struct bgen_vg *vg);
int get_ncombs(const struct bgen_vg *vg);

void close_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg);

void free(void *);

struct bgen_string string_duplicate(const struct bgen_string s);
