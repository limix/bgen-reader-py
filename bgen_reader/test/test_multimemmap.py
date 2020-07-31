import io
import pytest
from bgen_reader._environment import BGEN_READER_CACHE_HOME
from bgen_reader._multimemmap import MultiMemMap


def example_filepath3(filename):
    filepath = BGEN_READER_CACHE_HOME / "test_data" / filename
    if filepath.exists():
        filepath.unlink()
    return filepath


def test_errors():  # !!!cmk be sure these are run
    with pytest.raises(ValueError):
        MultiMemMap(
            example_filepath3("errors.mmm"), mode="w"
        )  # No support for mode 'w' ('w+' is OK)


def test_reads():
    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(write_file, mode="w+")
    assert (
        len(mmm_wplus) == 0
    )
    mm0 = mmm_wplus.append_empty("mm0", shape=(2, 3), dtype="int32")
    assert len(mmm_wplus) == 1
    mm0[:, :] = [[0, 1, 2], [3, 4, 5]]#!!!cmk test that with 'r' mode we can't change values
    mmm_wplus.append_empty("mm1", shape=(0, 3), dtype="float32")
    assert len(mmm_wplus) == 2
    assert mm0[0, 0] == 0
    del mmm_wplus

    with pytest.raises(FileNotFoundError):
        MultiMemMap(example_filepath3("doesntexisit.mmm"), mode="r")

    other_file = example_filepath3("otherfile.txt")
    with other_file.open('w') as fp:
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
    mmm_rplus.append_empty("mm3", shape=(1, 3), dtype="str")
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
    assert len(mmm_rplus)==2
    del mmm_r2
    del mmm_rplus


# !!!cmk test that 'F' and 'C' order are respected
# !!!cmk test appending name that is already there
# !!!cmk test a dtype of np.int32


def test_writes():
    write_file = example_filepath3("write.mmm")#!!!cmk test opening an existing file for write and having to re-started
    mmm_wplus = MultiMemMap(
        write_file, mode="w+", wplus_memmap_max=2, wplus_memmap_param_dtype="<U1"
    )
    with pytest.raises(ValueError):
        mmm_wplus.append_empty("0", shape=(3), dtype="str")  # Type is too long
    del mmm_wplus

    mmm_wplus = MultiMemMap(
        write_file, mode="w+", wplus_memmap_max=2, wplus_memmap_param_dtype="<U1"
    )
    mmm_wplus.append_empty("0", shape=(3), dtype="S")
    mmm_wplus.append_empty("1", shape=(3), dtype="S")
    with pytest.raises(ValueError):
        mmm_wplus.append_empty("2", shape=(5), dtype="S")  # Too many memmaps
    del mmm_wplus

    mmm_wplus = MultiMemMap(
        write_file, mode="r+", wplus_memmap_max=200, wplus_memmap_param_dtype="<U100"
    )  # memmap params are ignored with r+ and r
    with pytest.raises(ValueError):
        mmm_wplus.append_empty("mm2", shape=(1, 3), dtype="S")  # Too many memmaps
    del mmm_wplus

    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(write_file, mode="w+")
    mmm_wplus.append_empty("mm0", shape=(1, 3), dtype="str")
    with pytest.raises(KeyError):
        mmm_wplus.append_empty("mm0", shape=(1, 3), dtype="str")  # Duplicate name
    del mmm_wplus

    write_file = example_filepath3("write.mmm")
    mmm_wplus = MultiMemMap(write_file, mode="w+")
    zero_length = write_file.stat().st_size
    mmm_wplus.append_empty("mm0", shape=(1, 3), dtype="str")
    one_length = write_file.stat().st_size
    mmm_wplus.append_empty("mm1", shape=(1, 3), dtype="str")
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
    if True:
        test_errors()
        test_reads()
        test_writes()

    pytest.main([__file__])
