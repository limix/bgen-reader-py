from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(r"""
    typedef unsigned char byte;
    typedef int_fast64_t  inti;

    typedef struct BGenFile BGenFile;
    BGenFile* reader_open(const char *filepath);
    inti      reader_close(BGenFile *bgenfile);
    inti      reader_nsamples(BGenFile *bgenfile);
    inti      reader_nvariants(BGenFile *bgenfile);
    inti      read_variants(BGenFile *bgenfile, byte **ids);
""")

ffibuilder.set_source(
    "bgen_reader._ffi",
    r"""
    #include "bgen_reader/bgen_reader.h"

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

    inti read_variants(BGenFile *bgenfile, byte **ids)
    {
        VariantIdBlock *head_ref;
        inti e = bgen_reader_read_variantid_blocks(bgenfile, &head_ref);

        if (e != EXIT_SUCCESS)
            return EXIT_FAILURE;

        inti nvariants = bgen_reader_nvariants(bgenfile);
        inti i;

        for (i = 0; i < nvariants; ++i)
        {
            ids[i] = bgen_reader_strndup(head_ref->id, head_ref->id_length);
        }

        return EXIT_SUCCESS;
    }
""",
    libraries=['bgen_reader'],
    library_dirs=[join(get_config_var('prefix'), 'lib')],
    include_dirs=[join(get_config_var('prefix'), 'include')])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
