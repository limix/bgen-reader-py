import numpy as np


class MultiMemMap:  # !!!should be record and offer order 'F' vs 'C'?

    _expected_magic_string = "MultiMemMapV0"
    _bootstrap_dtype = "<U50"
    _bootstrap_max = 8
    _memmap_param_max = 8

    def __init__(self, filename, mode, metameta_dtype="<U50",  memmap_max = 25,
):
        # !!!cmk check all values of mode
        self._filename = filename
        self._mode = mode
        self._name_to_memmap = {}
        if filename.exists():
            self._offset = 0  # !!!cmk how come not using the inputted mode here?
            self._bootstrap = np.memmap(
                filename,
                dtype=self._bootstrap_dtype,
                mode="r+",
                offset=self._offset,
                shape=(self._bootstrap_max),
            )
            self._offset += self._bootstrap.size * self._bootstrap.itemsize
            assert self._magic_string == self._expected_magic_string, "Invalid file format. (Didn't find expected magic string.)"
            assert self._memmap_count <= self._memmap_max, "real assert"

            self._memmap_param = np.memmap(
                filename,
                dtype=metameta_dtype,
                mode="r+",
                offset=self._offset,
                shape=(self._memmap_max, self._memmap_param_max),
            )
            self._offset += self._memmap_param.size * self._memmap_param.itemsize

            for memmap_index in range(self._memmap_count):
                memmap = np.memmap(
                    filename,
                    dtype=self._get_memmap_dtype(memmap_index),
                    mode="r+",
                    offset=self._offset,
                    shape=self._get_memmap_shape(memmap_index),
                )
                self._offset += memmap.size * memmap.itemsize
                self._name_to_memmap[self._get_memmap_name(memmap_index)] = memmap
        else:
            assert mode == "w+", "MultiMemMap doesn't exist and mode isn't 'w+'"
            self._offset = 0
            self._bootstrap = np.memmap(
                self._filename,
                dtype=self._bootstrap_dtype,
                mode="w+",
                offset=self._offset,
                shape=(self._bootstrap_max),
            )
            self._offset += self._bootstrap.size * self._bootstrap.itemsize
            self._magic_string = self._expected_magic_string
            self._memmap_count = 0
            self._memmap_max = memmap_max
            self._bootstrap.flush() #!!!cmk offer (and use a global flush)

            self._memmap_param = np.memmap(
                self._filename,
                dtype=metameta_dtype,
                mode="r+",
                offset=self._offset,
                shape=(self._memmap_max, self._memmap_param_max),
            )
            self._offset += self._memmap_param.size * self._memmap_param.itemsize

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

    def _get_memmap_name(self,index):
        assert 0 <= index and index < self._memmap_count > 0, "Expect index between 0 (inclusive) and memmap_count (exclusive)"
        return self._memmap_param[index, 0]

    def _set_memmap_name(self, index, value):
        assert 0 <= index and index < self._memmap_count > 0, "Expect index between 0 (inclusive) and memmap_count (exclusive)"
        self._memmap_param[index, 0] = value

    def _get_memmap_dtype(self, index):
        assert 0 <= index and index < self._memmap_count > 0, "Expect index between 0 (inclusive) and memmap_count (exclusive)"
        return self._memmap_param[index, 1]

    def _set_memmap_dtype(self, index, value):
        assert 0 <= index and index < self._memmap_count > 0, "Expect index between 0 (inclusive) and memmap_count (exclusive)"
        self._memmap_param[index, 1] = str(value)  # cmk repr???

    def _get_memmap_shape(self, index):
        assert 0 <= index and index < self._memmap_count > 0, "Expect index between 0 (inclusive) and memmap_count (exclusive)"
        shape_as_str = self._memmap_param[index, 2]
        return tuple([int(num_as_str) for num_as_str in shape_as_str.split(",")])

    def _set_memmap_shape(self, index, value):
        assert 0 <= index and index < self._memmap_count > 0, "Expect index between 0 (inclusive) and memmap_count (exclusive)"
        self._memmap_param[index, 2] = str(value) #cmk repre???

    def __len__(self):
        assert len(self._name_to_memmap)==self._memmap_count,"real assert"
        return self._memmap_count

    def __getitem__(self, name):
        return self._name_to_memmap[name]

    def append_empty(
        self, name, shape, dtype
    ):  # !!!cmk say that all these dtypes must be strings, not types
        assert self._mode == "w+", "Can only append with mode 'w+'"
        assert self._memmap_count+1 < self._memmap_max, "The MultiMemMap contains no room for an additional memmap."
        self._memmap_count += 1
        self._set_memmap_name(self._memmap_count-1, name)
        self._set_memmap_dtype(self._memmap_count-1, dtype)
        self._set_memmap_shape(self._memmap_count-1, shape)
        self._memmap_param.flush()
        memmap = np.memmap(
            self._filename, dtype=dtype, mode="r+", offset=self._offset, shape=shape
        )
        # !!!cmk raise an error if every go over slots[2] e.g., 20
        self._bootstrap.flush()
        self._name_to_memmap[name] = memmap
        self._offset += memmap.size * memmap.itemsize
        memmap.flush()
        return memmap

    def popitem(
        self,
    ):  # As of Python 3.7 popitem remove the last item from a dictionary
        assert self._mode == "w+", "Can only popitem with mode 'w+'"
        assert self._memmap_count > 0, "The MultiMemMap contains no items to pop."
        name = self._get_memmap_name(self._memmap_count-1)
        self._set_memmap_name(self._memmap_count-1, None)
        self._memmap_count += -1
        self._memmap_param.flush()
        self._bootstrap.flush()
        memmap = self._name_to_memmap.pop(name)
        self._offset -= memmap.size * memmap.itemsize
        memmap._mmap.close()

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        if (
            hasattr(self, "_name_to_memmap") and self._name_to_memmap is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            for memmap in self._name_to_memmap.values():
                memmap._mmap.close()
            self._name_to_memmap.clear()
            del (
                self._name_to_memmap
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect.
        if (
            hasattr(self, "_bootstrap") and self._bootstrap is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            self._bootstrap._mmap.close()
            del (
                self._bootstrap
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect.
        if (
            hasattr(self, "_memmap_param") and self._memmap_param is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            self._memmap_param._mmap.close()
            del (
                self._memmap_param
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect
        if self._filename.stat().st_size > self._offset:
            with open(self._filename, "a") as fp:
                fp.truncate(self._offset)

    def __del__(self):
        self.__exit__()


# !!!cmk needs tests
