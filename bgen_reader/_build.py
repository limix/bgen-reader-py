from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(r"""

    typedef struct BGenFile BGenFile;
    BGenFile* open(const char *filepath);
""")

ffibuilder.set_source(
    "bgen_reader._ffi",
    r"""
    #include "bgen_reader/bgen_reader.h"

    BGenFile* open(const char *filepath)
    {
        return bgen_reader_open(filepath);
    }
""",
    libraries=['bgen_reader'],
    library_dirs=[join(get_config_var('prefix'), 'lib')],
    include_dirs=[join(get_config_var('prefix'), 'include')])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
