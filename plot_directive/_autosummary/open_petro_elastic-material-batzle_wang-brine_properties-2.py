import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_xlim(20, 350)
ax.set_xlabel("TEMPERATURE (°C)")
ax.set_ylabel("BULK MODULUS (GPa)")
xs = np.linspace(0.0, 350, 1000)
xs_short = np.linspace(0.0, 100, 1000)
#
def plot_brine_bulk_modulus(salinity, temperatures, pressure, *args, **kwargs):
    return ax.plot(
        temperatures,
        [
            bw.brine(x, pressure * 1e6, salinity * 1e6).bulk_modulus / 1e9
            for x in temperatures
        ],
        *args,
        **kwargs
    )
#
plot_brine_bulk_modulus(
    0.0,
    xs_short,
    0.1,
    "k--",
    label="PPM = 0",
)
plot_brine_bulk_modulus(
    0.15,
    xs_short,
    0.1,
    "k-",
    label="PPM = 150000",
)
plot_brine_bulk_modulus(
    0.3,
    xs_short,
    0.1,
    "k-.",
    label="PPM = 300000",
)
#
plot_brine_bulk_modulus(0.0, xs, 50, "k--")
plot_brine_bulk_modulus(0.15, xs, 50, "k-")
plot_brine_bulk_modulus(0.3, xs, 50, "k-.")
plot_brine_bulk_modulus(0.0, xs, 100, "k--")
plot_brine_bulk_modulus(0.15, xs, 100, "k-")
plot_brine_bulk_modulus(0.3, xs, 100, "k-.")
#
ax.legend()
