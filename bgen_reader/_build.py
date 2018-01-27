import os
import platform
from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

folder = os.path.dirname(os.path.abspath(__file__))


def windows_include_dirs():
    include_dirs = []
    if 'INCLUDE' in os.environ:
        include_dirs += [os.environ['INCLUDE']]
    if 'LIBRARY_INC' in os.environ:
        include_dirs += [os.environ['LIBRARY_INC']]
    if 'ProgramW6432' in os.environ:
        fld = join(os.environ['ProgramW6432'], 'bgen', 'include')
        if os.path.exists(fld):
            include_dirs += [fld]
    if 'ProgramFiles' in os.environ:
        fld = join(os.environ['ProgramFiles'], 'bgen', 'include')
        if os.path.exists(fld):
            include_dirs += [fld]
    return include_dirs


with open(join(folder, 'interface.h'), 'r') as f:
    ffibuilder.cdef(f.read())

with open(join(folder, 'interface.c'), 'r') as f:
    libraries = ['bgen', 'z']
    include_dirs = [join(get_config_var('prefix'), 'include')]
    library_dirs = [join(get_config_var('prefix'), 'lib')]

    if platform.system() == 'Windows':
        libraries += ['libzstd']
        include_dirs += windows_include_dirs()
    else:
        libraries += ['zstd']

    # INCLUDE

    ffibuilder.set_source(
        "bgen_reader._ffi",
        f.read(),
        libraries=libraries,
        library_dirs=include_dirs,
        include_dirs=library_dirs,
        language='c')

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
