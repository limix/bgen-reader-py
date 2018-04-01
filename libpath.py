"""
libpath
=======

Common path resolution for library and header files.

For instance, in Windows, the environment variable ``INCLUDE``, if exists,
is used to retrieve directories for header inclusion.
In Unix systems, ``/usr/local/lib`` is often used to store libraries.
In Windows, a library named ``z`` might exist under the name of ``libz``.
This module is a helper for handling those issues.

Example
-------

    >>> import platform
    >>> from os.path import join
    >>> from libpath import Windows, Unix
    >>>
    >>> if platform.system() == 'Windows':
    ...     s = Windows()
    ...     f = s.get_programfiles()
    ...     s.add_library_dir(join(f, 'zstd', 'lib'))
    ...     s.add_include_dir(join(f, 'zstd', 'include'))
    ...     libs = [s.find_libname('zstd')]
    ... else:
    ...     s = Unix()
    ...     libs = ['zstd']
    >>>
    >>> include_dirs = s.get_include_dirs()
    >>> library_dirs = s.get_library_dirs()
"""
import os
import struct
from os.path import exists, join
from sysconfig import get_config_var

__version__ = "0.0.1"


def bits_arch():
    """Determines the number of bits of the Python process.

        Return ``32`` or ``64``.
    """
    return struct.calcsize("P") * 8


class System(object):
    def __init__(self):
        self._library_dirs = []
        self._include_dirs = []

    def add_library_dir(self, d):
        self._library_dirs.append(d)

    def add_include_dir(self, d):
        self._include_dirs.append(d)


class Windows(System):

    def get_programfiles(self):
        if bits_arch() == 32:
            return self.get_32bits_programfiles()
        return self.get_64bits_programfiles()

    def get_64bits_programfiles(self):
        n = 'PROGRAMW6432'
        if n not in os.environ:
            raise RuntimeError("{} variable is not defined.".format(n))

        f = os.environ[n]
        if not exists(f):
            raise RuntimeError("Could not find {}.".format(f))

        return f

    def get_32bits_programfiles(self):
        n = 'PROGRAMFILES'
        if n not in os.environ:
            raise RuntimeError("{} variable is not defined.".format(n))

        f = os.environ[n]
        if not exists(f):
            raise RuntimeError("Could not find {}.".format(f))

        return f

    def get_include_dirs(self):
        dirs = [join(get_config_var('prefix'), 'include')]

        names = ['INCLUDE', 'LIBRARY_INC']
        vals = [os.environ[n] for n in names if n in os.environ]
        dirs += [d for v in vals for d in v.split(';')]

        dirs = [d for d in dirs if len(d) > 0 and exists(d)]
        return self._include_dirs + dirs

    def get_library_dirs(self):
        dirs = [join(get_config_var('prefix'), 'lib')]

        names = ['LIBRARY_LIB']
        vals = [os.environ[n] for n in names if n in os.environ]
        dirs += [d for v in vals for d in v.split(';')]

        dirs = [d for d in dirs if len(d) > 0 and exists(d)]
        return self._library_dirs + dirs

    def find_libname(self, name):
        """Try to infer the correct library name."""
        names = ["{}.lib", "lib{}.lib", "{}lib.lib"]
        names = [n.format(name) for n in names]
        dirs = self.get_library_dirs()
        for d in dirs:
            for n in names:
                if exists(join(d, n)):
                    return n[:-4]
        msg = "Could not find the {} library.".format(name)
        raise ValueError(msg)

    def __str__(self):
        msg = "Arch: {}bits\n".format(bits_arch())
        if bits_arch() == 32:
            msg += "ProgramFiles: {}\n".format(self.get_32bits_programfiles())
        if bits_arch() == 64:
            msg += "ProgramFiles32: {}\n".format(
                self.get_32bits_programfiles())
            msg += "ProgramFiles64: {}\n".format(
                self.get_64bits_programfiles())
        msg += "Include dirs: {}\n".format(self.get_include_dirs())
        msg += "Library dirs: {}\n".format(self.get_library_dirs())
        return msg


class Unix(System):

    def get_include_dirs(self):
        dirs = [join(get_config_var('prefix'), 'include')]
        dirs += ['/usr/include', '/usr/local/include']
        dirs = [d for d in dirs if len(d) > 0 and exists(d)]
        return self._include_dirs + dirs

    def get_library_dirs(self):
        dirs = [join(get_config_var('prefix'), 'lib')]
        dirs += ['/usr/lib', '/usr/local/lib']
        dirs = [d for d in dirs if len(d) > 0 and exists(d)]
        return self._library_dirs + dirs

    def __str__(self):
        msg = "Arch: {}bits\n".format(bits_arch())
        msg += "Include dirs: {}\n".format(self.get_include_dirs())
        msg += "Library dirs: {}\n".format(self.get_library_dirs())
        return msg
