import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_xlabel("PRESSURE (°C)")
ax.set_ylabel("VELOCITY (km/s)")
ax.set_ylim(0.9, 1.7)
ax.set_xlim(0.0, 40.0)
xs = np.linspace(1.0, 40, 1000)
xs_short = np.linspace(5.0, 40, 1000)
#
def plot_dead_oil_primary_velocity(temperature, pressures, *args, **kwargs):
   return ax.plot(
       pressures,
       [
           bw.dead_oil(temperature, x * 1e6, 916).primary_velocity / 1e3
           for x in pressures
       ],
       *args,
       **kwargs
   )
#
def plot_oil_primary_velocity(temperature, pressures, gor, *args, **kwargs):
   return ax.plot(
       pressures,
       [
           bw.oil(temperature, x * 1e6, 916, gor, 0.8).primary_velocity / 1e3
           for x in pressures
       ],
       *args,
       **kwargs
   )
#
plot_dead_oil_primary_velocity(
    22.8,
    xs,
    "k--",
)
plot_dead_oil_primary_velocity(
    72.0,
    xs,
    "k--",
)
#
plot_oil_primary_velocity(
   22.8,
   xs_short,
   85.0,
   "k--",
)
plot_oil_primary_velocity(
   72.0,
   xs_short,
   85.0,
   "k--",
)
#
ax.legend()
