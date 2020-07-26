import numpy as np

class MultiMemMap:
    def __init__(self, filename, mode, slot_dtype = 'int32', slots_shape =(3), metameta_dtype = '<U50'):
        #!!!cmk check all values of mode
        self._filename = filename
        self._mode = mode
        self._name_to_memmap = {}
        if filename.exists():
            self._offset = 0
            self._slots = np.memmap(filename,dtype=slot_dtype,mode='r+',offset=self._offset,shape=slots_shape)
            self._offset += self._slots.size * self._slots.itemsize
            slot_used,slot_count,self._metameta_count = self._slots
            assert slot_used<=slot_count

            self._metameta = np.memmap(filename,dtype=metameta_dtype,mode='r+',offset=self._offset,shape=(self._slots[1],self._slots[2]))
            self._offset += self._metameta.size * self._metameta.itemsize

            names = self._metameta[:slot_used,0]
            dtypes = self._metameta[:slot_used,1]
            shapes = []
            for shape_as_str in self._metameta[:slot_used,2]:
                shapes.append(tuple([int(num_as_str) for num_as_str in shape_as_str.split(',')]))

            for slot_index in range(slot_used):
                memmap = np.memmap(filename,dtype=dtypes[slot_index],mode='r+',offset=self._offset,shape=shapes[slot_index])
                self._offset += memmap.size * memmap.itemsize
                self._name_to_memmap[names[slot_index]]=memmap
        else:
            assert mode=="w+", "MultiMemMap doesn't exist and mode isn't 'w+'"
            self._offset = 0
            self._slots = np.memmap(self._filename,dtype=slot_dtype,mode='w+',offset=self._offset,shape=slots_shape)
            self._offset += self._slots.size * self._slots.itemsize
            self._slots[0] = 0
            self._slots[1] = 20
            self._slots[2] = 3
            self._slots.flush()

            self._metameta = np.memmap(self._filename,dtype=metameta_dtype,mode='r+',offset=self._offset,shape=(self._slots[1],self._slots[2]))
            self._offset += self._metameta.size * self._metameta.itemsize


    def __len__(self):
        return len(self._name_to_memmap)

    def __getitem__(self,name):
        return self._name_to_memmap[name]

    def append(self,name,nparray):
        assert self._mode == 'w+', "Can only append with mode 'w+'"
        slot_index = len(self)
        self._metameta[slot_index,0] = name
        self._metameta[slot_index,1] = str(nparray.dtype) #cmk repr???
        self._metameta[slot_index,2] = ','.join([str(i) for i in nparray.shape])
        self._metameta.flush()
        memmap = np.memmap(self._filename,dtype=nparray.dtype,mode='r+',offset=self._offset,shape=nparray.shape)
        self._offset += memmap.size * memmap.itemsize
        np.copyto(dst=memmap,src=nparray)
        memmap.flush()
        self._slots[0] = slot_index+1
        self._slots.flush()
        self._name_to_memmap[name] = memmap

    def append_copier(self, name, shape, dtype, copier):
        memmap = self.append_empty(name, shape, dtype)
        copier(memmap)
        memmap.flush()

    def append_empty(self, name, shape, dtype):
        assert self._mode == 'w+', "Can only append with mode 'w+'"
        slot_index = len(self)
        self._metameta[slot_index,0] = name
        self._metameta[slot_index,1] = str(dtype) #cmk repr???
        self._metameta[slot_index,2] = str(shape)
        self._metameta.flush()
        memmap = np.memmap(self._filename,dtype=dtype,mode='r+',offset=self._offset,shape=shape)
        self._slots[0] = slot_index+1
        self._slots.flush()
        self._name_to_memmap[name] = memmap
        self._offset += memmap.size * memmap.itemsize
        memmap.flush()
        return memmap


    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        if (
            hasattr(self, "_name_to_memmap")
            and self._name_to_memmap is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            for name, memmap in self._name_to_memmap.items():
                memmap.flush()
            self._name_to_memmap.clear()
            del (
                self._name_to_memmap
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect.
        if (
            hasattr(self, "_slots")
            and self._slots is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            self._slots.flush()
            del (
                self._slots
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect.
        if (
            hasattr(self, "_metameta")
            and self._metameta is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            self._metameta.flush()
            del (
                self._metameta
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect
    def __del__(self):
        self.__exit__()


#!!!cmk needs tests