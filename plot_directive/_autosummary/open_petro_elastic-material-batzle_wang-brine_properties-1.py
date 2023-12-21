import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_ylim(0.7, 1.3)
ax.set_xlim(20, 351)
ax.set_xlabel("TEMPERATURE (°C)")
ax.set_ylabel("DENSITY (g/cc)")
xs = np.linspace(0.0, 350, 1000)
#
def plot_brine_density(salinity, pressure, *args, **kwargs):
    return ax.plot(
        xs,
        [bw.brine(x, pressure * 1e6, salinity).density / 1000 for x in xs],
        *args,
        **kwargs
    )
#
plot_brine_density(0.0, 9.81, "k--")
plot_brine_density(0.0, 49, "k--")
plot_brine_density(0.0, 98.1, "k--")
#
for salinity in [0.02e6, 0.15e6, 0.24e6]:
    plot_brine_density(salinity, 9.81, "k-")
    plot_brine_density(salinity, 49, "k-")
    plot_brine_density(salinity, 98.1, "k-")
