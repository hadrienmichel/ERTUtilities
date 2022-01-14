# import pygimli as pg
from pygimli.physics.ert import ERTManager #, createGeometricFactors

ert = ERTManager(data='./VSL_Data/Line1_DDN6.dat', verbose=True)
mod = ert.invert()
ert.showResultAndFit()
ert.showModel(mod)
