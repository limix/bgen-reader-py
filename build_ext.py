import os
import platform
from os.path import join

from cffi import FFI

from libpath import Unix, Windows

ffibuilder = FFI()
libs = ['bgen', 'z', 'zstd', 'athr']

folder = os.path.dirname(os.path.abspath(__file__))

with open(join(folder, 'bgen_reader', 'interface.h'), 'r') as f:
    ffibuilder.cdef(f.read())

if platform.system() == 'Windows':
    s = Windows()
    f = s.get_programfiles()
    for lib in libs:
        s.add_library_dir(join(f, lib, 'lib'))
        s.add_include_dir(join(f, lib, 'include'))

    libs = [s.find_libname(lib) for lib in libs]
else:
    s = Unix()

ffibuilder.set_source(
    "bgen_reader._ffi",
    "#include \"bgen.h\"",
    libraries=libs,
    library_dirs=s.get_library_dirs(),
    include_dirs=s.get_include_dirs(),
    language='c')

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
