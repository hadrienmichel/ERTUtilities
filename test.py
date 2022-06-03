# %% EMI calibration:
import numpy as np
from matplotlib import pyplot as plt
dY = 0.1
spacings = [0.32, 0.71, 1.18]
Y_lim = np.arange(0, -5, -dY)
yCenters = (Y_lim[:-1]+Y_lim[1:])/2
### Horizontal loops calibration:
def sensitivityH(depth, spacing):
    return 2 - (4*depth/spacing)/(np.sqrt(4*(depth/spacing)**2 + 1))

calibrationH = np.zeros_like(spacings)
fig, ax = plt.subplots(1,1)
for j, space in enumerate(spacings):
    ax.plot([sensitivityH(-d, space) for d in yCenters], yCenters)
    # for i, d in enumerate(yCenters):
    #     calibrationH[j] += dY * res[i] * sensitivityH(d, space)

ax.set_xlabel('Sensitivity [/]')
ax.set_ylabel('Depth [m]')

print()

### Vertical loops calibration:
def sensitivityV(depth, spacing):
    return (4*depth/spacing)/((4*(depth/spacing)**2 + 1)**(3/2))

calibrationV = np.zeros_like(spacings)
fig, ax = plt.subplots(1,1)
for j, space in enumerate(spacings):
    ax.plot([sensitivityV(-d, space) for d in yCenters], yCenters)
    # for i, d in enumerate(yCenters):
    #     calibrationV[j] += dY * res[i] * sensitivityV(d, space)

ax.set_xlabel('Sensitivity [/]')
ax.set_ylabel('Depth [m]')

plt.show()