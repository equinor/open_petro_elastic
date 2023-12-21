import matplotlib.pyplot as plt
import numpy as np
#
import open_petro_elastic.material.batzle_wang as bw
fig, ax = plt.subplots()
ax.set_ylim(1.1, 1.8)
ax.set_xlim(1.1, 0.69)
ax.set_xlabel("ρ₀")
ax.set_ylabel("VELOCITY (km/sec)")
xs = np.linspace(1.1, 0.7, 1000)
ax.plot(
    xs,
    [bw.dead_oil(15.1, 0.1e6, x * 1e3).primary_velocity / 1000 for x in xs],
    "k-",
)
