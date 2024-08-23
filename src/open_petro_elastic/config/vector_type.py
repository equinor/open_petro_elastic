import numpy as np


class _ArrayMeta(type):
    def __getitem__(self, t):
        return type("Array", (Array,), {"__dtype__": t})


class Array(np.ndarray, metaclass=_ArrayMeta):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate_type

    @classmethod
    def validate_type(cls, val):
        assert val is not None
        dtype = getattr(cls, "__dtype__", None)
        if isinstance(dtype, tuple):
            dtype, shape = dtype
        else:
            shape = ()

        result = np.array(val, dtype=dtype, copy=False, ndmin=len(shape))
        assert not shape or len(shape) == len(result.shape)  # ndmin guarantees this

        if any(
            (shape[i] != -1 and shape[i] != result.shape[i]) for i in range(len(shape))
        ):
            result = result.reshape(shape)
        return result
