import io
from pathlib import Path
from typing import Tuple, Union

import numpy as np


class MultiMemMap:

    _expected_magic_string = "MultiMemMapV0"
    _bootstrap_dtype = "<U50"
    _bootstrap_max = 8
    _memmap_param_max = 8
    _mode_list = {"r", "r+", "w+"}
    # LATER document the meaning of each of these

    def __init__(
        self,
        filename: Union[str, Path],
        mode: str,
        wplus_memmap_param_dtype: str = "<U50",
        wplus_memmap_max: int = 25,  # LATER document the last two params are only used when mode is 'w+'
    ):
        self._filename = Path(filename)

        if mode not in self._mode_list:
            raise ValueError(f"invalid mode: '{mode}'")
        self._mode = mode
        if mode == "w+":
            self._create_new(wplus_memmap_param_dtype, wplus_memmap_max)
        else:
            self._read_existing()

    def _create_new(self, wplus_memmap_param_dtype: str, wplus_memmap_max: int):
        assert self._mode == "w+", "real assert"

        self._offset = 0
        self._bootstrap = np.memmap(
            self._filename,
            dtype=self._bootstrap_dtype,
            mode="w+",  # create file
            offset=self._offset,
            shape=(self._bootstrap_max),
            order="C",
        )
        self._offset += self._bootstrap.size * self._bootstrap.itemsize

        self._magic_string = self._expected_magic_string
        self._memmap_count = 0
        self._memmap_max = wplus_memmap_max
        self._memmap_param_dtype = wplus_memmap_param_dtype

        self._memmap_param = np.memmap(
            self._filename,
            dtype=self._memmap_param_dtype,
            mode="r+",  # file now exists, so 'r+' not 'w+'
            offset=self._offset,
            shape=(self._memmap_max, self._memmap_param_max),
            order="C",
        )
        self._offset += self._memmap_param.size * self._memmap_param.itemsize

        self._name_to_memmap = {}

    def _read_existing(self):
        if not self._filename.exists():
            raise FileNotFoundError(f"No such file or directory: '{self._filename}'")

        self._offset = 0
        self._bootstrap = np.memmap(
            self._filename,
            dtype=self._bootstrap_dtype,
            mode=self._mode,
            offset=self._offset,
            shape=(self._bootstrap_max),
            order="C",
        )
        self._offset += self._bootstrap.size * self._bootstrap.itemsize

        if self._magic_string != self._expected_magic_string:
            raise Exception("Invalid file format. (Didn't find expected magic string.)")
        assert self._memmap_count <= self._memmap_max, "real assert"

        self._memmap_param = np.memmap(
            self._filename,
            dtype=self._memmap_param_dtype,
            mode=self._mode,
            offset=self._offset,
            shape=(self._memmap_max, self._memmap_param_max),
            order="C",
        )
        self._offset += self._memmap_param.size * self._memmap_param.itemsize

        self._name_to_memmap = {}
        for memmap_index in range(self._memmap_count):
            memmap = np.memmap(
                self._filename,
                dtype=self._get_memmap_dtype(memmap_index),
                mode=self._mode,
                offset=self._offset,
                shape=self._get_memmap_shape(memmap_index),
                order=self._get_memmap_order(memmap_index),
            )
            self._offset += memmap.size * memmap.itemsize
            self._name_to_memmap[self._get_memmap_name(memmap_index)] = memmap

    def flush(self):
        self._bootstrap.flush()
        self._memmap_param.flush()
        for memmap in self._name_to_memmap.values():
            memmap.flush()

    @property
    def _magic_string(self):
        return self._bootstrap[0]

    @_magic_string.setter
    def _magic_string(self, value):
        self._bootstrap[0] = value

    @property
    def _memmap_count(self):
        return int(self._bootstrap[1])

    @_memmap_count.setter
    def _memmap_count(self, value):
        self._bootstrap[1] = str(value)

    @property
    def _memmap_max(self):
        return int(self._bootstrap[2])

    @_memmap_max.setter
    def _memmap_max(self, value):
        self._bootstrap[2] = str(value)

    @property
    def _memmap_param_dtype(self):
        return self._bootstrap[3]

    @_memmap_param_dtype.setter
    def _memmap_param_dtype(self, value):
        self._bootstrap[3] = value

    def _check_index(self, index):
        assert 0 <= index and index < self._memmap_count, "real assert"

    def _get_memmap_name(self, index):
        self._check_index(index)
        return self._memmap_param[index, 0]

    def _set_memmap_name(self, index, value):
        self._check_index(index)
        self._memmap_param[index, 0] = value
        if self._memmap_param[index, 0] != value:
            raise ValueError(f"Cannot save value as {self._memmap_param_dtype}")

    def _get_memmap_dtype(self, index):
        self._check_index(index)
        return self._memmap_param[index, 1]

    def _set_memmap_dtype(self, index, value):
        self._check_index(index)
        str_value = str(value)
        self._memmap_param[index, 1] = str_value
        if self._memmap_param[index, 1] != value:
            raise ValueError(f"Cannot save value as {self._memmap_param_dtype}")

    def _get_memmap_shape(self, index):
        self._check_index(index)
        shape_as_str = self._memmap_param[index, 2]
        return tuple((int(num_as_str) for num_as_str in shape_as_str.split(",")))

    def _set_memmap_shape(self, index, value):
        self._check_index(index)
        try:
            str_value = ",".join((str(num) for num in value))
        except TypeError:
            str_value = str(value)
        self._memmap_param[index, 2] = str_value
        if self._memmap_param[index, 2] != str_value:
            raise ValueError(f"Cannot save value as {self._memmap_param_dtype}")

    def _get_memmap_order(self, index):
        self._check_index(index)
        return self._memmap_param[index, 3]

    def _set_memmap_order(self, index, value):
        self._check_index(index)
        self._memmap_param[index, 3] = value
        if self._memmap_param[index, 3] != value:
            raise ValueError(f"Cannot save value as {self._memmap_param_dtype}")

    def __len__(self) -> int:
        assert len(self._name_to_memmap) == self._memmap_count, "real assert"
        return self._memmap_count

    def __contains__(self, item) -> bool:
        return item in self._name_to_memmap

    def __getitem__(self, name: str) -> np.memmap:
        return self._name_to_memmap[name]

    def append_empty(
        self,
        name: str,
        shape: Tuple[int],
        dtype: str,
        order: str = "C",
    ) -> np.memmap:  # Document that these dtypes must be strings, not types
        if self._mode not in {"r+", "w+"}:
            raise io.UnsupportedOperation("not writable")
        if self._memmap_count >= self._memmap_max:
            raise ValueError(
                "The MultiMemMap contains no room for an additional memmap."
            )
        if name in self._name_to_memmap:
            raise KeyError(f"A memmap with name '{name}' already exists")
        if order not in {"F", "C"}:
            raise TypeError("order not understood")

        self._memmap_count += 1
        self._set_memmap_name(self._memmap_count - 1, name)
        self._set_memmap_dtype(self._memmap_count - 1, dtype)
        self._set_memmap_shape(self._memmap_count - 1, shape)
        self._set_memmap_order(self._memmap_count - 1, order)

        memmap = np.memmap(
            self._filename,
            dtype=dtype,
            mode="r+",  # Because the file already exists
            offset=self._offset,
            shape=shape,
            order=order,
        )
        self._offset += memmap.size * memmap.itemsize
        self._name_to_memmap[name] = memmap

        if memmap.itemsize == 0:
            self.popitem()
            raise ValueError(
                f"Cannot use dtype '{dtype}' because it has a variable- or zero-size element"
            )

        return memmap

    def popitem(
        self,
    ):  # As of Python 3.7 popitem removes the last item from a dictionary
        if self._mode not in {"r+", "w+"}:
            raise io.UnsupportedOperation("not writable")
        if self._memmap_count == 0:
            raise KeyError("poptiem(): MultiMemMap is empty")

        name = self._get_memmap_name(self._memmap_count - 1)
        # Not necessary, but mark 'name' out (leave other params as is)
        self._set_memmap_name(self._memmap_count - 1, "")

        self._memmap_count -= 1
        memmap = self._name_to_memmap.pop(name)
        memmap._mmap.close()
        self._offset -= memmap.size * memmap.itemsize

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        # Python doesn't guarantee that __init__ was fully run, so check attrs.
        # We also remove attr. That allows __del__ and __exit__ to be called twice with no bad effect
        if hasattr(self, "_name_to_memmap") and self._name_to_memmap is not None:
            for memmap in self._name_to_memmap.values():
                memmap._mmap.close()
            self._name_to_memmap.clear()
            del self._name_to_memmap
        if hasattr(self, "_bootstrap") and self._bootstrap is not None:
            self._bootstrap._mmap.close()
            del self._bootstrap
        if hasattr(self, "_memmap_param") and self._memmap_param is not None:
            self._memmap_param._mmap.close()
            del self._memmap_param

        # If the file is longer than needed (because of 'popitem'), shorten it.
        if (
            hasattr(self, "_mode")
            and self._mode in {"r+", "w+"}
            and hasattr(self, "_filename")
            and self._filename.exists()
            and hasattr(self, "_offset")
            and self._filename.stat().st_size > self._offset
        ):
            with open(self._filename, "a") as fp:
                fp.truncate(self._offset)

    def __del__(self):
        self.__exit__()
