import os
import platform
from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

folder = os.path.dirname(os.path.abspath(__file__))

with open(join(folder, 'interface.h'), 'r') as f:
    ffibuilder.cdef(f.read())

with open(join(folder, 'interface.c'), 'r') as f:
    libraries = ['bgen', 'z']
    if platform.system() == 'Windows':
        libraries.append('libzstd')
    else:
        libraries.append('zstd')
    ffibuilder.set_source(
        "bgen_reader._ffi",
        f.read(),
        libraries=libraries,
        library_dirs=[join(get_config_var('prefix'), 'lib')],
        include_dirs=[join(get_config_var('prefix'), 'include')],
        language='c')

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
