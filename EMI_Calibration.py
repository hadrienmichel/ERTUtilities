from pygimli.physics import ert# import pygimli as pg
from pygimli import meshtools as mt
import pygimli as pg
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

mgr = ert.ERTManager(data='./TestFiles/Calibration.dat', verbose=True)
data = mgr.data
dX = 0.125
dY = 0.1
data['k'] = ert.createGeometricFactors(data)
data['rhoa'] = data['r'] * data['k']
mod = mgr.invert(data = data, paraMaxCellSize=0.015)
mgr.showResultAndFit()
mgr.showModel(mod)

m = pg.Mesh(mgr.paraDomain)
m['Res'] = mgr.paraModel(mod)
m.exportVTK('tmp.vtk')

### Extract column at the center of the profil:

dY = 0.1
x = 4
dXDevice = 1.5
vtkMesh = pv.read('tmp.vtk')
centers = vtkMesh.cell_centers()
val = centers.get_array('Res')
pts = centers.points
Y_up = vtkMesh.bounds[3]
Y_down = vtkMesh.bounds[2]
Y_lim = np.arange(0, -5, -dY)
yCenters = (Y_lim[:-1]+Y_lim[1:])/2
res = np.zeros_like(yCenters)
for j, lim in enumerate(yCenters):
    nbCells = 0
    for i, ptsCurr in enumerate(pts):
        if np.abs(ptsCurr[1]-lim) < dY/2 and np.abs(ptsCurr[0]-x) < dXDevice/2:
            nbCells += 1
            res[j] += val[i]
    res[j] /= nbCells

fig, ax = plt.subplots(1,1)
ax.plot(1000*np.divide(1, res), yCenters)
ax.set_xlabel('Conductivity [mS/m]')
ax.set_ylabel('Depth [m]')


# %% EMI calibration:
import numpy as np
from matplotlib import pyplot as plt
dY = 0.1
spacings = [0.32, 0.71, 1.18]
Y_lim = np.arange(0, -5, -dY)
yCenters = np.divide((Y_lim[:-1]+Y_lim[1:]),2)
### Horizontal loops calibration:
def sensitivityH(depth, spacing):
    return 2 - (4*depth/spacing)/(np.sqrt(4*(depth/spacing)**2 + 1))

calibrationH = np.zeros_like(spacings)
fig, ax = plt.subplots(1,1)
for j, space in enumerate(spacings):
    sens = [sensitivityH(-d, space) for d in yCenters]
    ax.plot(sens, yCenters)
    integ = np.sum(np.asarray(sens)*dY)
    for i, d in enumerate(yCenters):
        calibrationH[j] += dY * res[i] * sensitivityH(-d, space) / integ

ax.set_xlabel('Sensitivity [/]')
ax.set_ylabel('Depth [m]')

print(1000*1/calibrationH)

### Vertical loops calibration:
def sensitivityV(depth, spacing):
    return (4*depth/spacing)/((4*(depth/spacing)**2 + 1)**(3/2))

calibrationV = np.zeros_like(spacings)
fig, ax = plt.subplots(1,1)
for j, space in enumerate(spacings):
    sens = [sensitivityV(-d, space) for d in yCenters]
    ax.plot(sens, yCenters)
    integ = np.sum(np.asarray(sens)*dY)
    for i, d in enumerate(yCenters):
        calibrationV[j] += dY * res[i] * sensitivityV(-d, space) / integ

ax.set_xlabel('Sensitivity [/]')
ax.set_ylabel('Depth [m]')

print(1000*1/calibrationV)

plt.show()
