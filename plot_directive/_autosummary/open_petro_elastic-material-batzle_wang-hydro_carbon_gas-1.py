import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_ylim(0.0, 0.6)
#
ax.set_xlabel("TEMPERATURE (°C)")
ax.set_ylabel("GAS DENSITY (g/cm³)")
xs = np.linspace(0, 350, 1000)
#
def plot_gas_density(gas_gravity, pressure, *args, **kwargs):
    return ax.plot(
        xs,
        [bw.gas(x, pressure * 1e6, gas_gravity).density / 1e3 for x in xs],
        *args,
        **kwargs
    )
#
plot_gas_density(0.6, 10, "k--", label="G = 0.6")
plot_gas_density(1.2, 10, "k-", label="G = 1.2")
plot_gas_density(1.2, 25, "k-")
plot_gas_density(1.2, 50, "k-")
plot_gas_density(0.6, 25, "k--")
plot_gas_density(0.6, 50, "k--")
#
ax.text(15.0, 0.09, "10", rotation=-20)
ax.text(15.0, 0.21, "25", rotation=-32)
ax.text(15.0, 0.265, "50", rotation=-21)
#
ax.text(105.0, 0.18, "10", rotation=-31)
ax.text(100.0, 0.35, "25", rotation=-32)
ax.text(100.0, 0.415, "50MPa", rotation=-21)
#
ax.legend()
