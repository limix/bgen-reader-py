import pytest
from pathlib import Path


def test_rst_examples():
    import doctest
    for path in (Path(__file__).parent / "../../docs").glob("*.rst"):
        doctest.testfile(str(path))


if __name__ == "__main__":  #!!!cmk99 remove?
    pytest.main([__file__])
