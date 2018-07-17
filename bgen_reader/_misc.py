import errno
import stat
import sys
import os
from os.path import exists
from ._ffi import ffi

PY3 = sys.version_info >= (3,)

if not PY3:
    FileNotFoundError = IOError


def _group_readable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IRGRP)


def make_sure_bytes(p):
    try:
        p = p.encode()
    except AttributeError:
        pass
    return p


def make_sure_str(p):
    try:
        p = p.decode()
    except AttributeError:
        pass
    return p


def create_string(v):
    s = ffi.new("char[]", v.len)
    ffi.memmove(s, v.str, v.len)
    return ffi.string(s, v.len).decode()


def check_file_exist(filepath):
    if not exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)


def check_file_readable(filepath):
    if not _group_readable(filepath):
        msg = "You don't have file"
        msg += " permission for reading {}.".format(filepath)
        raise RuntimeError(msg)
