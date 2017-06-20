import os
from os.path import join
from sysconfig import get_config_var

import numpy
from cffi import FFI

ffibuilder = FFI()

folder = os.path.dirname(os.path.abspath(__file__))

with open(join(folder, 'interface.h'), 'r') as f:
    ffibuilder.cdef(f.read())

with open(join(folder, 'interface.c'), 'r') as f:
    ffibuilder.set_source(
        "bgen_reader._ffi",
        f.read(),
        libraries=['bgen', 'numpy'],
        library_dirs=[join(get_config_var('prefix'), 'lib')],
        include_dirs=[
            join(get_config_var('prefix'), 'include'),
            numpy.get_include()
        ],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
