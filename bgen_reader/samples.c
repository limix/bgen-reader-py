#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef MAX
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#endif

static void read_samples_part1(struct bgen_samples* samples, uint32_t nsamples,
                               uint32_t* samples_max_len)
{
    *samples_max_len = 0;

    for (uint32_t i = 0; i < nsamples; ++i) {
        struct bgen_string const* sample = bgen_samples_get(samples, i);
        *samples_max_len = MAX(*samples_max_len, (uint32_t)bgen_string_length(sample));
    }
}

static void read_samples_part2(struct bgen_samples const* samples, uint32_t nsamples,
                               char* const sample_array, uint32_t samples_stride)
{
    for (uint32_t i = 0; i < nsamples; ++i) {
        struct bgen_string const* sample = bgen_samples_get(samples, i);
        memcpy(sample_array + i * samples_stride, bgen_string_data(sample),
               bgen_string_length(sample));
    }
}
