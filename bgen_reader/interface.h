struct bgen_string;
struct bgen_string const* bgen_string_create(char const* data, size_t length);
void                      bgen_string_destroy(struct bgen_string const* bgen_string);
char const*               bgen_string_data(struct bgen_string const* bgen_string);
size_t                    bgen_string_length(struct bgen_string const* bgen_string);
bool                      bgen_string_equal(struct bgen_string a, struct bgen_string b);

struct bgen_file;
struct bgen_file*     bgen_file_open(char const* filepath);
void                  bgen_file_close(struct bgen_file const* bgen_file);
int32_t               bgen_file_nsamples(struct bgen_file const* bgen_file);
int32_t               bgen_file_nvariants(struct bgen_file const* bgen_file);
bool                  bgen_file_contain_samples(struct bgen_file const* bgen_file);
struct bgen_samples*  bgen_file_read_samples(struct bgen_file* bgen_file);
struct bgen_genotype* bgen_file_open_genotype(struct bgen_file* bgen_file,
                                              uint64_t          genotype_offset);

struct bgen_genotype;
void     bgen_genotype_close(struct bgen_genotype const* genotype);
int      bgen_genotype_read(struct bgen_genotype* genotype, double* probabilities);
uint16_t bgen_genotype_nalleles(struct bgen_genotype const* genotype);
bool     bgen_genotype_missing(struct bgen_genotype const* genotype, uint32_t index);
uint8_t  bgen_genotype_ploidy(struct bgen_genotype const* genotype, uint32_t index);
uint8_t  bgen_genotype_min_ploidy(struct bgen_genotype const* genotype);
uint8_t  bgen_genotype_max_ploidy(struct bgen_genotype const* genotype);
unsigned bgen_genotype_ncombs(struct bgen_genotype const* genotype);
bool     bgen_genotype_phased(struct bgen_genotype const* genotype);

struct bgen_metafile;
struct bgen_metafile* bgen_metafile_create(struct bgen_file* bgen_file, char const* filepath,
                                           uint32_t npartitions, int verbose);
struct bgen_metafile* bgen_metafile_open(char const* filepath);
uint32_t              bgen_metafile_npartitions(struct bgen_metafile const* metafile);
uint32_t              bgen_metafile_nvariants(struct bgen_metafile const* metafile);
struct bgen_partition const* bgen_metafile_read_partition(struct bgen_metafile const* metafile,
                                                          uint32_t partition);
int                          bgen_metafile_close(struct bgen_metafile const* metafile);

struct bgen_partition;
void                       bgen_partition_destroy(struct bgen_partition const* partition);
struct bgen_variant const* bgen_partition_get_variant(struct bgen_partition const* partition,
                                                      uint32_t                     index);
uint32_t                   bgen_partition_nvariants(struct bgen_partition const* partition);

struct bgen_samples;
void                      bgen_samples_destroy(struct bgen_samples const* samples);
struct bgen_string const* bgen_samples_get(struct bgen_samples const* samples, uint32_t index);

struct bgen_variant
{
    uint64_t                   genotype_offset;
    struct bgen_string const*  id;
    struct bgen_string const*  rsid;
    struct bgen_string const*  chrom;
    uint32_t                   position;
    uint16_t                   nalleles;
    struct bgen_string const** allele_ids;
};
