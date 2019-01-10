from ._ffi import ffi


def make_sure_bytes(p):
    try:
        p = p.encode()
    except AttributeError:
        pass
    return p


def create_string(v):
    s = ffi.new("char[]", v.len)
    ffi.memmove(s, v.str, v.len)
    return ffi.string(s, v.len).decode()


def bgen_str_to_str(s):
    if s.str == ffi.NULL:
        return ""
    return ffi.string(s.str, s.len).decode()
