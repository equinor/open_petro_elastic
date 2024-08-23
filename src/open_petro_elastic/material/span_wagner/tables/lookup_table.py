import functools

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def generate_lookup_table(func, x, y, filename):
    grids2d = np.meshgrid(x, y, indexing="ij")
    grids = [g.flatten() for g in grids2d]
    z = func(*grids)
    z_grid = z.reshape(grids2d[0].shape)
    np.savez(filename, z_grid=z_grid, x=x, y=y)
    return grids2d, z_grid


@functools.lru_cache(maxsize=10)
def load_lookup_table_interpolator(filename):
    data = np.load(filename)
    reg = RegularGridInterpolator(
        (data["x"], data["y"]), data["z_grid"], method="nearest", bounds_error=False
    )

    def _interp(_x, _y):
        _x, _y = np.asarray(_x), np.asarray(_y)
        pts = np.full((max(_x.size, _y.size), 2), fill_value=np.nan)
        pts[:, 0] = _x
        pts[:, 1] = _y
        res = reg(pts)
        if _x.ndim == 0 and _y.ndim == 0:
            return res[0]
        else:
            return res

    return _interp
