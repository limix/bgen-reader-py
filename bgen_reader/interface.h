struct bgen_str {
    int len;
    char *str;
};

/* Variant metadata. */
struct bgen_vm {
    long vaddr; /* variant offset-address */
    struct bgen_str id;
    struct bgen_str rsid;
    struct bgen_str chrom;
    int position; /* base-pair position */
    int nalleles; /* number of alleles */
    struct bgen_str *allele_ids;
};

struct bgen_file;
struct bgen_vi;
struct bgen_vg;

struct bgen_file *bgen_open(const char *);
void bgen_close(struct bgen_file *);

int bgen_nsamples(const struct bgen_file *);
int bgen_nvariants(const struct bgen_file *);
int bgen_contain_samples(const struct bgen_file *);

struct bgen_str *bgen_read_samples(struct bgen_file *, int);
void bgen_free_samples(const struct bgen_file *, struct bgen_str *);

struct bgen_mf *bgen_create_metafile(struct bgen_file *, const char *, int, int);
struct bgen_mf *bgen_open_metafile(const char *);
int bgen_metafile_nparts(struct bgen_mf *);
int bgen_metafile_nvars(struct bgen_mf *);
struct bgen_vm *bgen_read_partition(struct bgen_mf *, int, int *);
void bgen_free_partition(struct bgen_vm *, int);
int bgen_close_metafile(struct bgen_mf *);
int bgen_ncombs(const struct bgen_vg *);
int bgen_phased(const struct bgen_vg *);
int bgen_missing(const struct bgen_vg *, size_t);
int bgen_ploidy(const struct bgen_vg *, size_t);
struct bgen_vg *bgen_open_genotype(struct bgen_file *, long);
void bgen_close_genotype(struct bgen_vg *vg);
int bgen_read_genotype(struct bgen_file *, struct bgen_vg *, double *);

/* Deprecated functions. */
int bgen_create_variants_metadata_file(const char *, const char *, int);
struct bgen_var *bgen_read_variants_metadata(struct bgen_file *, struct bgen_vi **, int);
int bgen_store_variants_metadata(const struct bgen_file *, struct bgen_var *,
                                 struct bgen_vi *, const char *);
struct bgen_var *bgen_load_variants_metadata(const struct bgen_file *, const char *,
                                             struct bgen_vi **, int);
void bgen_close_variant_genotype(struct bgen_vi *, struct bgen_vg *);
void bgen_free_variants_metadata(const struct bgen_file *, struct bgen_var *);
void bgen_free_index(struct bgen_vi *);
struct bgen_vg *bgen_open_variant_genotype(struct bgen_vi *, size_t);
int bgen_read_variant_genotype(struct bgen_vi *, struct bgen_vg *, double *);
int bgen_sample_ids_presence(const struct bgen_file *);
int bgen_max_nalleles(struct bgen_vi *);
