import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
#
fig, ax = plt.subplots()
ax.set_ylim(0.5, 2.0)
ax.set_xlim(0, 350)
ax.set_xlabel("TEMPERATURE (°C)")
ax.set_ylabel("VELOCITY (km/s)")
xs = np.linspace(0.0, 350, 1000)
#
ax.plot(
    xs,
    [bw.water(x, 1.0).primary_velocity / 1000 for x in xs],
    "k-",
    xs,
    [bw.water(x, 50).primary_velocity / 1000 for x in xs],
    "k-",
    xs,
    [bw.water(x, 100).primary_velocity / 1000 for x in xs],
    "k-",
)
