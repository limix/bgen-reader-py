from __future__ import unicode_literals

import os

from bgen_reader import print_something


def test_bgen_reader():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    print_something(os.path.join(folder, b"example.32bits.bgen"))
