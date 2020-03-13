from ._ffi import ffi, lib


def make_sure_bytes(p):
    try:
        p = p.encode()
    except AttributeError:
        pass
    return p


def create_string(bgen_string):
    length = lib.bgen_string_length(bgen_string)
    c_string = ffi.new("char[]", length)
    ffi.memmove(c_string, lib.bgen_string_data(bgen_string), length)
    return ffi.string(c_string, length).decode()


def bgen_str_to_str(s):
    if s.str == ffi.NULL:
        return ""
    return ffi.string(s.str, s.len).decode()
