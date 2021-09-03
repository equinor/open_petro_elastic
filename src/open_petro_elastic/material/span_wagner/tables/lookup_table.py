import numpy as np
from scipy.interpolate import RectBivariateSpline


def generate_lookup_table(func, x, y, filename):
    grids2d = np.meshgrid(x, y, indexing='ij')
    grids = [g.flatten() for g in grids2d]
    z = func(*grids)
    z_grid = z.reshape(grids2d[0].shape)
    np.savez(filename, z_grid=z_grid, x=x, y=y)
    return grids2d, z_grid


def load_lookup_table_interpolator(filename):
    data = np.load(filename)
    rbs = RectBivariateSpline(data['x'], data['y'], data['z_grid'])
    box = np.min(data['x']), np.max(data['x']), np.min(data['y']), np.max(data['y'])

    def _interp(_x, _y):
        outside = (_x < box[0])
        outside |= (_x > box[1])
        outside |= (_y < box[2])
        outside |= (_y > box[3])
        res = np.full_like(outside, dtype=np.float, fill_value=np.nan)
        _x = _x if np.asarray(_x).size == 1 else _x[~outside]
        _y = _y if np.asarray(_y).size == 1 else _y[~outside]
        res[~outside] = rbs.ev(_x, _y)
        return res

    return _interp
