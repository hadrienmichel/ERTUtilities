'''Complete list of options for Inversion objects: 
(from: https://www.pygimli.org/gimliapi/classGIMLI_1_1Inversion.html)
    - setRelativeError(pg.Vector or float): set the relative error for all data points
                                            (float) or points-by-points (pg.Vector). 
                                            If set to 0.0, no data weigthing or correction.
                                            DEFAULT = 0.01 (1%)
    - setAbsoluteError(pg.Vector or float): set absolute error for all data points (float) 
                                            or points-by-points (pg.Vector).
                                            DEFAULT = not used
    - setTransData(transformation): set a transformation for tha data (log for example).
                                    DEFAULT = Identity
    - setTransModel(transformation): set a transformation for tha model (log for example).
                                     DEFAULT = Identity
    - setVerbose(bool): set the verbose status.
                        DEFAULT = False
    - setDoSave(bool): set debug output.
                       DEFAULT = False
    - setLineSearch(bool): set line search optimization of the update.
                           DEFAULT = True
    - setBlockyModel(bool): set the constraint on the model to L1 instead of L2.
                            DEFAULT = False
    - setRobustData(bool): set the constraint on the data to L1 instead of L2.
                           DEFAULT = False
    - setLambda(float): set the value of the regularization parameter.
                        DEFAULT = 50.0
    - setLambdaFactor(float): set the multiplier for the regularization at each iteration 
                              (see Marquardt scheme).
                              DEFAULT = 1.0
    - setLambdaMinimum(float): set minimum value of the regularization parameter.
                               DEFAULT = 1.0
    - setLocalRegularization(bool): set local regularization (true) or global (false).
                                    DEFAULT = FALSE
    - setOptimizeLambda(bool): optimize Lambda using the L-curve criteria.
                               DEFAULT = False
    - setMaxIter(int): set the maximum number of iterations for the inversion.
                       DEFAULT = 20
    - setMaxCGLSIter(int): set the maximum number of iterations for the conjugated
                           gradient least square optimization at each iteration.
                           DEFAULT = 200
    - setCGLSTolerance(float): set the tolerance for CGLS (default scaling = -1).
                               DEFAULT = -1.0
    - stopAtChi1(bool): set the inversion to stop at Chi²=1 (stopping criteria).
                        DEFAULT = True
    - setDeltaPhiAbortPercent(float): set the minimum acceptable delta in the objective 
                                      function between iterations (stopping criteria).
                                      DEFAULT = 2.0
    - setMarquardtScheme(float): set the inversion to follow a Marquardt scheme. The
                                 input value is the damping factor that multiplies 
                                 lambda at each iteration.
                                 DEFAULT = 0.8
    - setRecalcJacobian(bool): set the inversion to recalculate the Jacobian matrix
                               at each iteration.
                               DEFAULT = True
    - setBroydenUpdate(bool): set the update to Broyden (scaling by transform function
                              derivatives).
                              DEFAULT = False
    - setModel(pg.Vector): set the starting model for the inversion.
                           (default to mean value)
    - setReferenceModel(pg.Vector): set the reference model for the inversion
                                    (default to starting model).
    - setCWeight(pg.Vector): set the constraint weight vector (boundary control).
    - setMWeight(pg.Vector): set the model weigth for each cells.

This list of options is complete and adress directly the GIMLi API.

The list of transfoms is available at https://www.pygimli.org/gimliapi/classGIMLI_1_1Trans.html

'''

import pybert as pb
import pygimli as pg
import numpy as np


tdip=pb.TDIPdata(filename='AgatheData\PP9_ChessionAmont_topo.ohm', verbose=True)
data = tdip.data
mgr = pb.ERTManager(data = tdip.data)


#tdip.showRhoa()
mesh = pg.meshtools.readTetgen('D:\Documents\dox@uliege\Article HydroGPhy\Profiles GPHY\DataFiles\Inversion\mesh\mesh.1')

# mesh = mgr.createMesh(data, quality = 33.6,  paraDX = 0.5) # , paraDepth = 50, paraMaxCellSize = 1,


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

# tdip.showResistivity(cMap='jet')



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
# tdip.showChargeability(cMap='jet',cMin = 0, cMax = 25)


### Export VTK 
tdip.pd.addData("Res", tdip.res)
tdip.pd.addData("Coverage", tdip.coverage)
tdip.pd.addData("IP", tdip.m*1000)
tdip.pd.exportVTK('AgatheData\PP9_ChessionAmont')