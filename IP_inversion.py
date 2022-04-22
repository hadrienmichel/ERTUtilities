# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:57:24 2022

@author: tomde
"""


import pybert as pb
import pygimli as pg
import numpy as np


tdip=pb.TDIPdata(filename='AgatheData\PP9_ChessionAmont.ohm', verbose=True)
data = tdip.data
mgr = pb.ERTManager(data = tdip.data)


#tdip.showRhoa()


mesh = mgr.createMesh(data, quality = 33.6,  paraDX = 0.5) # , paraDepth = 50, paraMaxCellSize = 1,


# # Invert resistivity
lam = 20
verbose = True
robustData = False
chi1Opt = False # Optimize lambda for chi2 = 1
zWeight = 1.0
blockyModel = False

tdip.invertRhoa(mesh=mesh) # ,lam = 10, robustData= False, verbose = True)
ertMgr = tdip.ERT
inv = ertMgr.fw.inv
inv.setRobustData(robustData)
inv.setVerbose(verbose)
inv.setLambda(lam)
inv.setDeltaPhiAbortPercent(1)

tdip.res = inv.runChi1()

tdip.coverage = ertMgr.coverage()
tdip.pd = ertMgr.paraDomain

tdip.showResistivity(cMap='jet')



# # # Invert chargeability
ma = tdip.data['ip']*0.001
maerr = tdip.data['ipErr']*0.001
fop = tdip.ERT.fop
mesh = tdip.ERT.mesh
res = tdip.res
fIp = pb.tdip.mipmodelling.DCIPMModelling(fop, mesh, res)

if fop.regionManager().regionCount() > 1:
    fIp.region(1).setBackground(True)
    fIp.region(2).setConstraintType(1)

fIp.regionManager().setZWeight(zWeight)

INV = pg.core.RInversion(ma, fIp, True, False)
tD, tM = pg.core.RTrans(), pg.core.RTransLogLU(0.0001, 0.999)
INV.setTransData(tD)
INV.setTransModel(tM)
mstart = pg.Vector(len(res), np.abs(pg.median(ma)))
INV.setModel(mstart)
INV.setAbsoluteError(maerr)
INV.setLambda(lam)
INV.setRobustData(robustData)
INV.setBlockyModel(blockyModel)
INV.setDeltaPhiAbortPercent(1)

tdip.m = INV.run()
tdip.mafwd = INV.response()

# tdip.invertMa(nr = 0, lam= 10, ma = tdip.data['ip']*0.001)
tdip.showChargeability(cMap='jet',cMin = 0, cMax = 25)


### Export VTK 
tdip.pd.addData("Res", tdip.res)
tdip.pd.addData("Coverage", tdip.coverage)
tdip.pd.addData("IP", tdip.m*1000)
tdip.pd.exportVTK('AgatheData\PP9_ChessionAmont')