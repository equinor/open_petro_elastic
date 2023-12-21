import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_xlabel("TEMPERATURE (°C)")
ax.set_ylabel("OIL BULK MODULUS (MPa)")
ax.set_ylim(0.0, 3400)
xs = np.linspace(0.0, 350, 1000)
xs_short = np.linspace(0.0, 100, 1000)
#
def plot_oil_bulk_modulus(
    reference_density, temperatures, pressure, *args, **kwargs
):
    return ax.plot(
        temperatures,
        [
            bw.dead_oil(x, pressure * 1e6, reference_density * 1000).bulk_modulus
            / 1e6
            for x in temperatures
        ],
        *args,
        **kwargs
    )
#
plot_oil_bulk_modulus(
    1.00,
    xs_short,
    0.1,
    "k-.",
    label="10 deg. API",
)
plot_oil_bulk_modulus(
    0.88,
    xs_short,
    0.1,
    "k-",
    label="30 deg. API",
)
plot_oil_bulk_modulus(
    0.78,
    xs_short,
    0.1,
    "k--",
    label="50 deg. API",
)
plot_oil_bulk_modulus(0.78, xs, 25, "k--")
plot_oil_bulk_modulus(0.78, xs, 50, "k--")
plot_oil_bulk_modulus(0.88, xs, 25, "k-")
plot_oil_bulk_modulus(0.88, xs, 50, "k-")
plot_oil_bulk_modulus(1.0, xs, 25, "k-.")
plot_oil_bulk_modulus(1.0, xs, 50, "k-.")
#
ax.legend()
