# import os
# import platform
# from os.path import join

# from cffi import FFI

# from libpath import Unix, Windows

# ffibuilder = FFI()
# libs = ["bgen", "z", "zstd", "athr"]

# folder = os.path.dirname(os.path.abspath(__file__))

# with open(join(folder, "bgen_reader", "interface.h"), "r") as f:
#     ffibuilder.cdef(f.read())

# if platform.system() == "Windows":
#     s = Windows()
#     f = s.get_programfiles()
#     for lib in libs:
#         s.add_library_dir(join(f, lib, "lib"))
#         s.add_include_dir(join(f, lib, "include"))

#     libs = [s.find_libname(lib) for lib in libs]
# else:
#     s = Unix()

# library_dirs = s.get_library_dirs()
# extra_link_args = []
# if platform.system() == "Darwin":
#     if len(library_dirs) > 0:
#         extra_link_args += ["-Wl,-rpath," + ",-rpath,".join(library_dirs)]

# ffibuilder.set_source(
#     "bgen_reader._ffi",
#     '#include "bgen.h"',
#     libraries=libs,
#     library_dirs=library_dirs,
#     include_dirs=s.get_include_dirs(),
#     extra_link_args=extra_link_args,
#     language="c",
# )

# if __name__ == "__main__":
#     ffibuilder.compile(verbose=True)

import os
import platform
from os.path import join
from typing import List

from cffi import FFI

from libpath import System, Unix, Windows

ffibuilder = FFI()
libs = ["bgen", "z", "zstd", "athr"]

folder = os.path.dirname(os.path.abspath(__file__))

with open(join(folder, "bgen_reader", "interface.h"), "r") as f:
    ffibuilder.cdef(f.read())

if platform.system() == "Windows":
    win = Windows()
    progfiles = win.get_programfiles()
    for lib in libs:
        win.add_library_dir(join(progfiles, lib, "lib"))
        win.add_include_dir(join(progfiles, lib, "include"))

    libs = [win.find_libname(lib) for lib in libs]
    system: System = win
else:
    system = Unix()

library_dirs = system.get_library_dirs()
extra_link_args: List[str] = []
if platform.system() == "Darwin":
    if len(library_dirs) > 0:
        extra_link_args += ["-Wl,-rpath," + ",-rpath,".join(library_dirs)]

ffibuilder.set_source(
    "bgen_reader._ffi",
    r"""
    #include "bgen/bgen.h"
    """,
    libraries=libs,
    library_dirs=library_dirs,
    include_dirs=system.get_include_dirs(),
    extra_link_args=extra_link_args,
    language="c",
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
