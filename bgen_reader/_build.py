from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

with open('interface.h', 'r') as f:
    ffibuilder.cdef(f.read())

with open('interface.c', 'r') as f:
    ffibuilder.set_source(
        "bgen_reader._ffi",
        f.read(),
        libraries=['bgen'],
        library_dirs=[join(get_config_var('prefix'), 'lib')],
        include_dirs=[join(get_config_var('prefix'), 'include')])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
