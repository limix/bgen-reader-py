from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(r"""
    typedef int_fast64_t  inti;

    typedef struct BGenFile BGenFile;
    BGenFile* reader_open(const char *filepath);
    inti      reader_close(BGenFile *bgenfile);
    inti      reader_nsamples(BGenFile *bgenfile);
    inti      reader_nvariants(BGenFile *bgenfile);
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
""",
    libraries=['bgen_reader'],
    library_dirs=[join(get_config_var('prefix'), 'lib')],
    include_dirs=[join(get_config_var('prefix'), 'include')])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
