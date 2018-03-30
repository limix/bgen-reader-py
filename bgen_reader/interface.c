#include <stdio.h>
#include <stdlib.h>

#include "bgen.h"

struct bgen_file *open_bgen(const char *filepath) {
  return bgen_open(filepath);
}

void close_bgen(struct bgen_file *bgen) { bgen_close(bgen); }

int get_nsamples(const struct bgen_file *bgen) { return bgen_nsamples(bgen); }

int get_nvariants(const struct bgen_file *bgen) { return bgen_nvariants(bgen); }

int sample_ids_presence(const struct bgen_file *bgen) {
  return bgen_sample_ids_presence(bgen);
}

struct bgen_string *read_samples(struct bgen_file *bgen, int verbose) {
  return bgen_read_samples(bgen, verbose);
}

void free_samples(const struct bgen_file *bgen, struct bgen_string *samples) {
  bgen_free_samples(bgen, samples);
}

struct bgen_var *read_variants_metadata(struct bgen_file *bgen,
                                        struct bgen_vi **index, int verbose) {
  return bgen_read_variants_metadata(bgen, index, verbose);
}

void free_variants_metadata(const struct bgen_file *bgen,
                            struct bgen_var *variants) {
  bgen_free_variants_metadata(bgen, variants);
}

void free_index(struct bgen_vi *index) { bgen_free_index(index); }

struct bgen_vg *open_variant_genotype(struct bgen_vi *index,
                                      size_t variant_idx) {
  return bgen_open_variant_genotype(index, variant_idx);
}

void read_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg,
                           double *probabilities) {
  bgen_read_variant_genotype(index, vg, probabilities);
}

int get_nalleles(const struct bgen_vg *vg) { return bgen_nalleles(vg); }

int get_ploidy(const struct bgen_vg *vg) { return bgen_ploidy(vg); }

int get_ncombs(const struct bgen_vg *vg) { return bgen_ncombs(vg); }

void close_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg) {
  bgen_close_variant_genotype(index, vg);
}

struct bgen_string string_duplicate(struct bgen_string s) {
  struct bgen_string r;

  r.str = malloc(s.len);
  memcpy(r.str, s.str, s.len);
  r.len = s.len;
  return r;
}
