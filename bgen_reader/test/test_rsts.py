import pytest
from pathlib import Path


def test_rst_examples():
    import doctest

    for path in (Path(__file__).parent / "../../docs").glob("*.rst"):
        rel_path = path.relative_to(Path(__file__).parent)
        doctest.testfile(str(rel_path))


if __name__ == "__main__":
    pytest.main([__file__])
