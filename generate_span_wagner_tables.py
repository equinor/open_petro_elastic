"""
This script is used to generate the tables included in open_petro_elastic.material.span_wagner.tables. It serves as
documentation for how the tables were generated, as well as a template for generating tables of different resolution or
bounds.
"""
import os
import numpy as np
from open_petro_elastic.material.span_wagner.tables.lookup_table import generate_lookup_table
from open_petro_elastic.material.span_wagner.carbon_dioxide import carbon_dioxide_density
from time import perf_counter


t_density = 273.15 + np.hstack((np.arange(0, 50, step=0.05), np.arange(50, 200.1, step=0.5)))
p_density = np.hstack((np.arange(0.5, 10, step=0.01), np.arange(10, 100.1, step=0.1)))

out_dir = os.path.join('src', 'open_petro_elastic', 'material', 'span_wagner', 'tables')

t0 = perf_counter()
density = generate_lookup_table(
    lambda _t, _p: carbon_dioxide_density(_t, _p, interpolate=False, force_vapor='auto'),
    t_density, p_density, os.path.join(out_dir, 'carbon_dioxide_density.npz')
)
t1 = perf_counter()
print(f'Generated density table in {t1 - t0} seconds')
