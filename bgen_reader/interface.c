#include <stdio.h>
#include <stdlib.h>

#include "bgen.h"

struct bgen_file *open_bgen(const char *filepath) {
  return bgen_open(filepath);
}

void close_bgen(struct bgen_file *bgen) { bgen_close(bgen); }

int get_nsamples(struct bgen_file *bgen) { return bgen_nsamples(bgen); }

int get_nvariants(struct bgen_file *bgen) { return bgen_nvariants(bgen); }

bgen_string *read_samples(struct bgen_file *bgen) {
  return bgen_read_samples(bgen, 0);
}

void free_samples(const struct bgen_file *bgen, bgen_string *samples) {
  bgen_free_samples(bgen, samples);
}

struct bgen_var *read_variants(struct bgen_file *bgen, struct bgen_vi **index) {
  return bgen_read_variants_metadata(bgen, index, 0);
}

void free_variants(const struct bgen_file *bgen, struct bgen_var *variants) {
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

int get_nalleles(struct bgen_vg *vg) { return bgen_nalleles(vg); }

int get_ploidy(struct bgen_vg *vg) { return bgen_ploidy(vg); }

int get_ncombs(struct bgen_vg *vg) { return bgen_ncombs(vg); }

void close_variant_genotype(struct bgen_vi *index, struct bgen_vg *vg) {
  bgen_close_variant_genotype(index, vg);
}

bgen_string string_duplicate(bgen_string s) {
  bgen_string r;

  r.str = malloc(s.len);
  memcpy(r.str, s.str, s.len);
  r.len = s.len;
  return r;
}

int sample_ids_presence(struct bgen_file *bgen) {
  return bgen_sample_ids_presence(bgen);
}
