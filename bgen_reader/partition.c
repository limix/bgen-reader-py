#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef MAX
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#endif

static void read_partition_part1(struct bgen_partition const* partition, uint32_t* position,
                                 uint16_t* nalleles, uint64_t* offset, uint32_t* id_max_len,
                                 uint32_t* rsid_max_len, uint32_t* chrom_max_len,
                                 uint32_t* allele_ids_max_len)
{
    uint32_t nvariants = bgen_partition_nvariants(partition);

    *id_max_len = 0;
    *rsid_max_len = 0;
    *chrom_max_len = 0;
    *allele_ids_max_len = 0;

    for (uint32_t i = 0; i < nvariants; ++i) {
        struct bgen_variant const* v = bgen_partition_get_variant(partition, i);
        position[i] = v->position;
        nalleles[i] = v->nalleles;
        offset[i] = v->genotype_offset;

        *id_max_len = MAX(*id_max_len, (uint32_t)bgen_string_length(v->id));
        *rsid_max_len = MAX(*rsid_max_len, (uint32_t)bgen_string_length(v->rsid));
        *chrom_max_len = MAX(*chrom_max_len, (uint32_t)bgen_string_length(v->chrom));

        if (v->nalleles == 0)
            continue;

        uint32_t length = v->nalleles - 1;
        for (uint16_t j = 0; j < v->nalleles; ++j)
            length += (uint32_t)bgen_string_length(v->allele_ids[j]);
        *allele_ids_max_len = MAX(*allele_ids_max_len, length);
    }
}

static void read_partition_part2(struct bgen_partition const* partition, char* const id,
                                 uint32_t id_stride, char* const rsid, uint32_t rsid_stride,
                                 char* const chrom, uint32_t chrom_stride,
                                 char* const allele_ids, uint32_t allele_ids_stride)
{
    uint32_t nvariants = bgen_partition_nvariants(partition);
    for (uint32_t i = 0; i < nvariants; ++i) {
        struct bgen_variant const* v = bgen_partition_get_variant(partition, i);

        memcpy(id + i * id_stride, bgen_string_data(v->id), bgen_string_length(v->id));

        memcpy(rsid + i * rsid_stride, bgen_string_data(v->rsid), bgen_string_length(v->rsid));

        memcpy(chrom + i * chrom_stride, bgen_string_data(v->chrom),
               bgen_string_length(v->chrom));

        size_t j = 0;
        for (uint16_t r = 0; r < v->nalleles; ++r) {

            memcpy(allele_ids + i * allele_ids_stride + j, bgen_string_data(v->allele_ids[r]),
                   bgen_string_length(v->allele_ids[r]));
            j += bgen_string_length(v->allele_ids[r]);

            if (r + 1 < v->nalleles) {
                memcpy(allele_ids + i * allele_ids_stride + j, ",", 1);
                j += 1;
            }
        }
    }
}
