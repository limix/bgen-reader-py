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

struct bgen_file;
struct bgen_vi;
struct bgen_vg;

struct bgen_file *bgen_open(const char *filepath);
void bgen_close(struct bgen_file *bgen);

int bgen_nsamples(const struct bgen_file *bgen);
int bgen_nvariants(const struct bgen_file *bgen);
int bgen_sample_ids_presence(const struct bgen_file *bgen);

struct bgen_string *bgen_read_samples(struct bgen_file *bgen, int verbose);
void bgen_free_samples(const struct bgen_file *bgen,
                       struct bgen_string *samples);

struct bgen_var *bgen_read_variants_metadata(struct bgen_file *bgen,
                                             struct bgen_vi **index,
                                             int verbose);
void bgen_free_variants_metadata(const struct bgen_file *bgen,
                                 struct bgen_var *variants);

struct bgen_vg *bgen_open_variant_genotype(struct bgen_vi *index,
                                           size_t variant_idx);
void bgen_read_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg,
                                double *probabilities);
int bgen_nalleles(const struct bgen_vg *vg);
int bgen_missing(const struct bgen_vg *vg, size_t index);
int bgen_ploidy(const struct bgen_vg *vg, size_t index);
int bgen_min_ploidy(const struct bgen_vg *vg);
int bgen_max_ploidy(const struct bgen_vg *vg);
int bgen_ncombs(const struct bgen_vg *vg);
int bgen_phased(const struct bgen_vg *vg);
void bgen_close_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg);

int bgen_store_variants_metadata(const struct bgen_file *bgen,
                                 struct bgen_var *variants, struct bgen_vi *vi,
                                 const char *filepath);
struct bgen_var *bgen_load_variants_metadata(const struct bgen_file *bgen,
                                             const char *filepath,
                                             struct bgen_vi **vi, int verbose);
int bgen_create_variants_metadata_file(const char *bgen_fp, const char *vi_fp,
                                       int verbose);
