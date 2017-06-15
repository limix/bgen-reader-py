from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(r"""
    typedef unsigned char byte;
    typedef int_fast64_t  inti;

    typedef struct string
    {
        inti len;
        byte *s;
    } string;

    typedef struct BGenFile BGenFile;
    BGenFile* reader_open(const char *filepath);
    inti      reader_close(BGenFile *bgenfile);
    inti      reader_nsamples(BGenFile *bgenfile);
    inti      reader_nvariants(BGenFile *bgenfile);
    inti      reader_read_variants(BGenFile *bgenfile, string **ids,
                                   string **rsids, string **chroms,
                                   inti *position, inti *nalleles);
    void      free(void *ptr);
""")

ffibuilder.set_source(
    "bgen_reader._ffi",
    r"""
    #include <stdlib.h>

    #include "bgen_reader/bgen_reader.h"

    typedef struct string
    {
        inti len;
        byte *s;
    } string;

    BGenFile* reader_open(const char *filepath)
    {
        return bgen_reader_open(filepath);
    }

    inti reader_close(BGenFile *bgenfile)
    {
        return bgen_reader_close(bgenfile);
    }

    inti reader_nsamples(BGenFile *bgenfile)
    {
        return bgen_reader_nsamples(bgenfile);
    }

    inti reader_nvariants(BGenFile *bgenfile)
    {
        return bgen_reader_nvariants(bgenfile);
    }

    inti reader_read_variants(BGenFile *bgenfile, string **ids,
                              string **rsids, string **chroms,
                              inti *position, inti *nalleles)
    {
        VariantIdBlock *head_ref = NULL;
        VariantIdBlock *ref = NULL;
        inti e = bgen_reader_read_variantid_blocks(bgenfile, &head_ref);

        if (e != EXIT_SUCCESS)
            return EXIT_FAILURE;

        inti nvariants = bgen_reader_nvariants(bgenfile);
        inti i;

        for (i = 0; i < nvariants; ++i)
        {
            ids[i] = malloc(sizeof(string));
            ids[i]->len = head_ref->id_length;
            ids[i]->s = bgen_reader_strndup(head_ref->id, head_ref->id_length);

            rsids[i] = malloc(sizeof(string));
            rsids[i]->len = head_ref->rsid_length;
            rsids[i]->s = bgen_reader_strndup(head_ref->rsid,
                                              head_ref->rsid_length);

            chroms[i] = malloc(sizeof(string));
            chroms[i]->len = head_ref->chrom_length;
            chroms[i]->s = bgen_reader_strndup(head_ref->chrom,
                                               head_ref->chrom_length);

            position[i] = head_ref->position;
            nalleles[i] = head_ref->nalleles;

            ref = head_ref;
            head_ref = head_ref->next;

            e = bgen_reader_free_variantid_block(ref);
        }

        return EXIT_SUCCESS;
    }
""",
    libraries=['bgen_reader'],
    library_dirs=[join(get_config_var('prefix'), 'lib')],
    include_dirs=[join(get_config_var('prefix'), 'include')])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
