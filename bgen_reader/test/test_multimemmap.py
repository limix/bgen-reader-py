import io

import pytest

from bgen_reader._environment import BGEN_READER_CACHE_HOME
from bgen_reader._multimemmap import MultiMemMap


def example_filepath3(filename):
    filepath = BGEN_READER_CACHE_HOME / "test_data" / filename
    if filepath.exists():
        filepath.unlink()
    return filepath


def test_errors():
    with pytest.raises(ValueError):
        MultiMemMap(
            example_filepath3("errors.mmm"), mode="w"
        )  # No support for mode 'w' ('w+' is OK)


def test_reads():
    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(write_file, mode="w+")
    assert len(mmm_wplus) == 0
    mm0 = mmm_wplus.append_empty("mm0", shape=(2, 3), dtype="int32")
    assert len(mmm_wplus) == 1
    mm0[:, :] = [
        [0, 1, 2],
        [3, 4, 5],
    ]
    mmm_wplus.append_empty("mm1", shape=(0, 3), dtype="float32")
    assert len(mmm_wplus) == 2
    assert mm0[0, 0] == 0
    del mmm_wplus

    with pytest.raises(FileNotFoundError):
        MultiMemMap(example_filepath3("doesntexisit.mmm"), mode="r")

    other_file = example_filepath3("otherfile.txt")
    with other_file.open("w") as fp:
        fp.write("this is not a mmm file")

    with pytest.raises(Exception):
        MultiMemMap(other_file, mode="r")

    with MultiMemMap(write_file, mode="r") as mmm_r1:
        assert len(mmm_r1) == 2
        with pytest.raises(io.UnsupportedOperation):
            mmm_r1.append_empty("mm3", shape=(1, 3), dtype="str")
        with pytest.raises(io.UnsupportedOperation):
            mmm_r1.popitem()
        with pytest.raises(ValueError):
            mmm_r1["mm0"][1, 2] = -50
        mmm_r2 = MultiMemMap(write_file, mode="r")
        assert len(mmm_r2) == 2

    mmm_rplus = MultiMemMap(write_file, mode="r+")
    mmm_rplus.append_empty("mm3", shape=(1, 3), dtype="<U10")
    assert len(mmm_rplus) == 3
    mm0 = mmm_rplus["mm0"]
    assert mm0[1, 2] == 5
    mm0[1, 2] = -5
    assert mm0[1, 2] == -5
    mmm_rplus.flush()
    mmm_r3 = MultiMemMap(write_file, mode="r")
    assert len(mmm_r3) == 3
    assert mmm_r3["mm0"][1, 2] == -5
    del mmm_r3
    mmm_rplus.popitem()
    assert len(mmm_rplus) == 2
    del mmm_r2
    del mmm_rplus


def test_more():
    write_file = example_filepath3("write.mmm")
    with MultiMemMap(write_file, mode="w+") as _:
        pass
    zero_length = write_file.stat().st_size
    with MultiMemMap(write_file, mode="r+") as mmm_rplus:
        with pytest.raises(ValueError):
            mmm_rplus.append_empty("mm0", shape=(3, 5), dtype="str")
        mmm_rplus.append_empty("mm0", shape=(3, 5), dtype="<U10")
    one_length = write_file.stat().st_size
    assert zero_length < one_length  # Expect file to grow
    # Opening the file in w+ will create it from scratch again
    with MultiMemMap(write_file, mode="w+") as _:
        pass
    assert write_file.stat().st_size == zero_length

    with MultiMemMap(write_file, mode="r+") as mmm_rplus:
        mm0 = mmm_rplus.append_empty("mm0", shape=(3, 5), dtype="<U10")
        assert mm0.flags["C_CONTIGUOUS"]
        with pytest.raises(TypeError):
            mmm_rplus.append_empty("mm1", shape=(3, 5), dtype="<U10", order="K")
        mm1 = mmm_rplus.append_empty("mm1", shape=(3, 5), dtype="<U10", order="F")
        assert mm1.flags["F_CONTIGUOUS"]
    with MultiMemMap(write_file, mode="r") as mmm_r:
        assert mmm_r["mm0"].flags["C_CONTIGUOUS"]
        assert mmm_r["mm1"].flags["F_CONTIGUOUS"]

    with MultiMemMap(write_file, mode="w+") as mmm_wplus:
        mmm_wplus.append_empty("mm0", shape=(3, 5), dtype="<U10")
        assert "mm0" in mmm_wplus
        mmm_wplus.append_empty("mm1", shape=(3, 5), dtype="<U10")
        mmm_wplus.append_empty("mm2", shape=(3, 5), dtype="<U10")
        mmm_wplus.popitem()
        mmm_wplus.popitem()
        mmm_wplus.close()
    assert write_file.stat().st_size == one_length

    with MultiMemMap(
        write_file, mode="w+", wplus_memmap_param_dtype="<U10"
    ) as mmm_wplus:
        with pytest.raises(ValueError):
            mmm_wplus.append_empty(
                "toolong_toolong_toolong", shape=(3, 5), dtype="<U10"
            )
    with MultiMemMap(
        write_file, mode="w+", wplus_memmap_param_dtype="<U10"
    ) as mmm_wplus:
        with pytest.raises(ValueError):
            mmm_wplus.append_empty(
                "mm0", shape=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), dtype="<U10"
            )


def test_writes():
    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(
        write_file, mode="w+", wplus_memmap_max=2, wplus_memmap_param_dtype="<U1"
    )
    with pytest.raises(ValueError):
        mmm_wplus.append_empty("0", shape=(3), dtype="<U10")  # Type is too long
    del mmm_wplus

    mmm_wplus = MultiMemMap(
        write_file, mode="w+", wplus_memmap_max=2, wplus_memmap_param_dtype="<U5"
    )
    mmm_wplus.append_empty("0", shape=(3), dtype="<U10")
    mmm_wplus.append_empty("1", shape=(3), dtype="<U10")
    with pytest.raises(ValueError):
        mmm_wplus.append_empty("2", shape=(5), dtype="<U10")  # Too many memmaps
    del mmm_wplus

    mmm_wplus = MultiMemMap(
        write_file, mode="r+", wplus_memmap_max=200, wplus_memmap_param_dtype="<U100"
    )  # memmap params are ignored with r+ and r
    with pytest.raises(ValueError):
        mmm_wplus.append_empty("mm2", shape=(1, 3), dtype="<U10")  # Too many memmaps
    del mmm_wplus

    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(write_file, mode="w+")
    mmm_wplus.append_empty("mm0", shape=(1, 3), dtype="<U10")
    with pytest.raises(KeyError):
        mmm_wplus.append_empty("mm0", shape=(1, 3), dtype="<U10")  # Duplicate name
    del mmm_wplus

    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(write_file, mode="w+")
    zero_length = write_file.stat().st_size
    mmm_wplus.append_empty("mm0", shape=(1, 3), dtype="<U10")
    one_length = write_file.stat().st_size
    mmm_wplus.append_empty("mm1", shape=(1, 3), dtype="<U10")
    two_length = write_file.stat().st_size
    del mmm_wplus
    assert two_length == write_file.stat().st_size
    mmm_rplus = MultiMemMap(write_file, mode="r+")
    mmm_rplus.popitem()
    assert len(mmm_rplus) == 1
    del mmm_rplus
    assert one_length == write_file.stat().st_size  # On close, files get trimmed
    mmm_rplus = MultiMemMap(write_file, mode="r+")
    mmm_rplus.popitem()
    assert len(mmm_rplus) == 0
    with pytest.raises(KeyError):
        mmm_rplus.popitem()  # Can't pop if no items
    del mmm_rplus
    assert zero_length == write_file.stat().st_size


if __name__ == "__main__":
    pytest.main([__file__])
