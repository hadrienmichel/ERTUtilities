import numpy as np
import matplotlib.pylab as plt

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert

###############################################################################
# Create a measurement scheme for 64 electrodes with a spacing of 1m
scheme = ert.createERTData(
    elecs=np.linspace(start=-31.5, stop=31.5, num=64),
    schemeName='dd')

###############################################################################
# Mesh generation
world = mt.createWorld(
    start=[-50, 0], end=[50, -15], layers=[-0.3, -2.0, -3.2, -9.80], worldMarker=True)

anomaly = mt.createCircle(
    pos=[0, -1], radius=0.15, marker=6
)

innerAnomaly = mt.createCircle(pos=[0,-1], radius=0.14, marker=7)

plc = world + anomaly + innerAnomaly

# local refinement of mesh near electrodes
for s in scheme.sensors():
    plc.createNode(s + [0.0, -0.2])

mesh_coarse = mt.createMesh(plc, quality=33)
# additional refinements
mesh = mesh_coarse.createH2()

# Create a P2-optimized mesh (quadratic shape functions)
mesh = mesh.createP2()

pg.show(plc, marker=True)
###############################################################################
# Prepare the model parameterization
# We have two markers here: 1: background 2: circle anomaly
rhomap = [
    [1, 100],
    [2, 40],
    [3, 250],
    [4, 50],
    [5, 400],
    [6, 1e-1],
    [7, 25]
]

# For visualization, map the rhomap into the actual mesh
rho = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)
fig, axes = plt.subplots(1, 1, figsize=(16 / 2.54, 16 / 2.54))
pg.show(mesh, data=rho, ax=axes, label=r"$\rho$ [$\Omega$m]")
fig.tight_layout()

###############################################################################
# Do the actual forward modeling
data = ert.simulate(
    mesh,
    res=rhomap,
    scheme=scheme,
    noiseAbs=0.0,
    noiseLevel=0.05,
)

###############################################################################
# Visualize the modeled data

# Please note the apparent negative (resistivity) phases!
fig, axes = plt.subplots(1, 1, figsize=(16 / 2.54, 16 / 2.54))
ert.showERTData(data, vals=data['rhoa'], ax=axes)

fig.tight_layout()

###############################################################################
# Inversion of the data
mgr = ert.ERTManager(data)
inv = mgr.invert()
mgr.showResultAndFit()
mgr.showResult()

plt.show()
