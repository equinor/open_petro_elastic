import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_ylim(0.0, 650.0)
ax.set_xlabel("TEMPERATURE (°C)")
ax.set_ylabel("GAS BULK MODULUS (MPa)")
xs = np.linspace(0.0, 350, 1000)
xs_short = np.linspace(40, 350, 1000)
#
def plot_gas_bulk_modulus(gas_gravity, temperatures, pressure, *args, **kwargs):
    return ax.plot(
        temperatures,
        [
            bw.gas(x, pressure * 1e6, gas_gravity).bulk_modulus / 1e6
            for x in temperatures
        ],
        *args,
        **kwargs
    )
#
plot_gas_bulk_modulus(0.6, xs, 10, "k--", label="G = 0.6")
plot_gas_bulk_modulus(1.2, xs, 10, "k-", label="G = 1.2")
plot_gas_bulk_modulus(1.2, xs, 25, "k-")
plot_gas_bulk_modulus(1.2, xs_short, 50, "k-")
plot_gas_bulk_modulus(0.6, xs, 25, "k--")
plot_gas_bulk_modulus(0.6, xs, 50, "k--")
#
ax.text(100.0, 25, "10", rotation=0)
ax.text(70.0, 60, "25", rotation=-2)
ax.text(45.0, 120, "50", rotation=-18)
#
ax.text(8.0, 22, "10", rotation=-31)
ax.text(22.0, 300, "25", rotation=-70)
ax.text(50.0, 500, "50MPa", rotation=-70)
#
ax.legend()
