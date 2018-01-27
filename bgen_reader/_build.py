import os
import platform
from os.path import join
from sysconfig import get_config_var

from cffi import FFI

ffibuilder = FFI()

folder = os.path.dirname(os.path.abspath(__file__))


def windows_dirs(prefix, lib):
    dirs = []
    if 'PROGRAMW6432' in os.environ:
        fld = join(os.environ['PROGRAMW6432'], lib, prefix)
        if os.path.exists(fld):
            dirs += [fld]
    if 'PROGRAMFILES' in os.environ:
        fld = join(os.environ['PROGRAMFILES'], lib, prefix)
        if os.path.exists(fld):
            dirs += [fld]
    return dirs


def windows_include_dirs():
    include_dirs = []
    if 'INCLUDE' in os.environ:
        include_dirs += [os.environ['INCLUDE']]
    if 'LIBRARY_INC' in os.environ:
        include_dirs += [os.environ['LIBRARY_INC']]
    include_dirs += windows_dirs('include', 'bgen')
    include_dirs += windows_dirs('include', 'zlib')
    include_dirs += windows_dirs('include', 'zstd')
    include_dirs += windows_dirs('include', 'athr')
    return include_dirs


def windows_library_dirs():
    library_dirs = []
    if 'LIBRARY_LIB' in os.environ:
        library_dirs += [os.environ['LIBRARY_LIB']]
    library_dirs += windows_dirs('lib', 'bgen')
    library_dirs += windows_dirs('lib', 'zlib')
    library_dirs += windows_dirs('lib', 'zstd')
    library_dirs += windows_dirs('lib', 'athr')
    return library_dirs


def windows_find_libname(lib, library_dirs):
    names = [
        "{}.lib".format(lib), "lib{}.lib".format(lib), "{}lib.lib".format(lib)
    ]
    folders = [f for ldir in library_dirs for f in ldir.split(';')]
    for f in folders:
        for n in names:
            if os.path.exists(join(f, n)):
                return n[:-4]

    raise RuntimeError("{} library not found.")


with open(join(folder, 'interface.h'), 'r') as f:
    ffibuilder.cdef(f.read())

with open(join(folder, 'interface.c'), 'r') as f:

    include_dirs = [join(get_config_var('prefix'), 'include')]
    library_dirs = [join(get_config_var('prefix'), 'lib')]

    if platform.system() == 'Windows':
        include_dirs += windows_include_dirs()
        library_dirs += windows_library_dirs()
        libraries = [
            windows_find_libname('zstd', library_dirs),
            windows_find_libname('z', library_dirs),
            windows_find_libname('athr', library_dirs),
            windows_find_libname('bgen', library_dirs)
        ]
    else:
        libraries = ['bgen', 'athr', 'z', 'zstd']

    ffibuilder.set_source(
        "bgen_reader._ffi",
        f.read(),
        libraries=libraries,
        library_dirs=library_dirs,
        include_dirs=include_dirs,
        language='c')

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
