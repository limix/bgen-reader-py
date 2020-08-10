#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

static void read_ploidy(struct bgen_genotype const* genotype, uint8_t* ploidy,
                        uint32_t nsamples)
{
    for (uint32_t i = 0; i < nsamples; ++i)
        ploidy[i] = bgen_genotype_ploidy(genotype, i);
}

static void read_missing(struct bgen_genotype const* genotype, bool* missing,
                         uint32_t nsamples)
{
    for (uint32_t i = 0; i < nsamples; ++i)
        missing[i] = bgen_genotype_missing(genotype, i);
}
