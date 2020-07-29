import numpy as np


class MultiMemMap:  # !!!should be record and offer order 'F' vs 'C'?

    _bootstrap_dtype = "int32"
    _bootstrap_length = 3

    def __init__(self, filename, mode, metameta_dtype="<U50"):
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
                shape=(self._bootstrap_length),
            )
            self._offset += self._bootstrap.size * self._bootstrap.itemsize
            assert self._memmap_count <= self._memmap_max, "real assert"

            self._memmap_param = np.memmap(
                filename,
                dtype=metameta_dtype,
                mode="r+",
                offset=self._offset,
                shape=(self._memmap_max, self._metameta_count),
            )
            self._offset += self._memmap_param.size * self._memmap_param.itemsize

            try: #!!!cmk
                names = self._memmap_param[: self._memmap_count, 0]
            except:
                print('!!!cmk')
            dtypes = self._memmap_param[: self._memmap_count, 1]
            shapes = []
            for shape_as_str in self._memmap_param[: self._memmap_count, 2]:
                shapes.append(
                    tuple([int(num_as_str) for num_as_str in shape_as_str.split(",")])
                )

            for memmap_index in range(self._memmap_count):
                memmap = np.memmap(
                    filename,
                    dtype=dtypes[memmap_index],
                    mode="r+",
                    offset=self._offset,
                    shape=shapes[memmap_index],
                )
                self._offset += memmap.size * memmap.itemsize
                self._name_to_memmap[names[memmap_index]] = memmap
        else:
            assert mode == "w+", "MultiMemMap doesn't exist and mode isn't 'w+'"
            self._offset = 0
            self._bootstrap = np.memmap(
                self._filename,
                dtype=self._bootstrap_dtype,
                mode="w+",
                offset=self._offset,
                shape=(self._bootstrap_length),
            )
            self._offset += self._bootstrap.size * self._bootstrap.itemsize
            self._memmap_count = 0
            self._memmap_max = 20
            self._metameta_count = 3
            self._bootstrap.flush() #!!!cmk offer (and use a global flush)

            self._memmap_param = np.memmap(
                self._filename,
                dtype=metameta_dtype,
                mode="r+",
                offset=self._offset,
                shape=(self._memmap_max, self._metameta_count),
            )
            self._offset += self._memmap_param.size * self._memmap_param.itemsize

    @property
    def _memmap_count(self):
        return self._bootstrap[0]

    @_memmap_count.setter
    def _memmap_count(self, value):
        self._bootstrap[0] = value

    @property
    def _memmap_max(self):
        return self._bootstrap[1]

    @_memmap_max.setter
    def _memmap_max(self, value):
        self._bootstrap[1] = value

    @property
    def _metameta_count(self):
        return self._bootstrap[2]

    @_metameta_count.setter
    def _metameta_count(self, value):
        self._bootstrap[2] = value

    @property
    def _last_memmap_name(self):
        assert self._memmap_count > 0, "With no memmaps, 'last' is not defined."
        return self._memmap_param[self._memmap_count-1, 0]

    @_last_memmap_name.setter
    def _last_memmap_name(self, value):
        assert self._memmap_count > 0, "With no memmaps, 'last' is not defined."
        self._memmap_param[self._memmap_count-1, 0] = value

    @property
    def _last_memmap_dtype(self):
        assert self._memmap_count > 0, "With no memmaps, 'last' is not defined."
        return self._memmap_param[self._memmap_count-1, 1]

    @_last_memmap_dtype.setter
    def _last_memmap_dtype(self, value):
        assert self._memmap_count > 0, "With no memmaps, 'last' is not defined."
        self._memmap_param[self._memmap_count-1, 1] = value

    @property
    def _last_memmap_shape(self):
        assert self._memmap_count > 0, "With no memmaps, 'last' is not defined."
        return self._memmap_param[self._memmap_count-1, 2]

    @_last_memmap_shape.setter
    def _last_memmap_shape(self, value):
        assert self._memmap_count > 0, "With no memmaps, 'last' is not defined."
        self._memmap_param[self._memmap_count-1, 2] = value

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
        self._last_memmap_name = name
        self._last_memmap_dtype = str(dtype)  # cmk repr???
        self._last_memmap_shape = str(shape)
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
        assert self._mode == "w+", "Can only append with mode 'w+'"
        name = self._last_memmap_name
        self._last_memmap_name = None
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
