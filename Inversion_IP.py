#!/usr/bin/env python
# coding: utf-8
# Author: David Caterina - ULiège (2022)
# # Inversion of IP data

# In[63]:


# Common imports
import numpy as np
from matplotlib import pyplot as plt
import pygimli as pg
from pygimli.physics.ert import createGeometricFactors, ERTManager
from pygimli.meshtools import createMesh
import pybert as pb
import os.path
from copy import deepcopy

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# ### Defining some usefull functions

# In[64]:


def plotCompareResults(data, invResults):
    relError = ((abs(data['rhoa']-invResults))/data['rhoa'] * 100) # Relative error in % for each data point
    figComp = plt.figure(figsize=(10,5))
    ax2 = figComp.add_subplot(122)
    ax2.plot(data['rhoa'], invResults, 'bo', label='Full Data')
    ax2.axline((0,0), slope=1, color='k', label='Unity line')
    ax2.set_xlabel('Measured Apparent resistivities [Ohm.m]')
    ax2.set_ylabel('Computed Apparent resistivities [Ohm.m]')
    ax1 = figComp.add_subplot(121)
    ax1.hist(relError, color='b')
    ax1.set_xlabel('Relative error [%]')
    ax1.set_ylabel('Number of data points')
    plt.plot()
    return


# In[65]:


def plotCompareResultsIP(data, invResults):
    ma=data['ip']*0.001
    relError = ((abs(ma-invResults))/ma * 100) # Relative error in % for each data point
    figComp = plt.figure(figsize=(10,5))
    ax2 = figComp.add_subplot(122)
    ax2.plot(ma, invResults, 'bo', label='Full Data')
    ax2.axline((0,0), slope=1, color='k', label='Unity line')
    ax2.set_xlabel('Measured Apparent chargeability [mV/V]')
    ax2.set_ylabel('Computed Apparent resistivities [mV/V]')
    ax1 = figComp.add_subplot(121)
    ax1.hist(relError, color='b')
    ax1.set_xlabel('Relative error [%]')
    ax1.set_ylabel('Number of data points')
    plt.plot()
    return


# ### File system and names

# In[66]:

plt.close('all')
directory = './P1/Pygimli/Orozco-slater/'
dataFile = 'P1.ohm'
saveDir = 'Results/'
if os.path.exists(directory+"TrimmedDataset.ohm"): # check if TrimmedDataset.ohm exists
    os.remove(directory+"TrimmedDataset.ohm")

if not os.path.isdir(directory+saveDir):
    os.mkdir(directory+saveDir)
    
xlabel="Distance (m)"
ylabel="Elevation (m)"


# ### Base inversion and parameters

# In[68]:


# Loading the dataset
tdip = pb.TDIPdata(directory+dataFile)
data = tdip.data

# meshInit = createMesh(data.sensors())
if not data.allNonZero('k'):
    data['k'] = createGeometricFactors(data, numerical=True)
if not data.allNonZero('rhoa'):
    data['rhoa'] = data['r']*data['k']

# Remove negative apparent resistivity
dataLength0=len(data["rhoa"])
data.remove(data["rhoa"]<=0)
dataLength1=len(data["rhoa"])

print(f'{dataLength0-dataLength1} data with negative resistivity, removed from the dataset')

if not data.allNonZero('err'):
    absoluteError = 0.001
    relativeError = 0.03
    data['err'] = relativeError + absoluteError / data['rhoa'] # In % of the resistance

dataErr=deepcopy(np.asarray((data['err'])))

    
mgr = pb.ERTManager(data)
# Mesh creation
paraDepth = 20 # Depth to obtain at any points on the profile
quality = 33.6
paraMaxCellSize = 0.5
paraDX = 0.3
# Changes to the paraDepth 
positions = data.sensors()
meshInit = mgr.createMesh(data, paraDepth=paraDepth, quality=quality-2, paraMaxCellSize=paraMaxCellSize*4, paraDX=paraDX*3)
mgr.setMesh(meshInit)
pg.show(mgr.paraDomain)


# ### Inversion of the data (round 1)

# In[69]:


# # Invert resistivity
lam = 20
verbose = True
robustData = True
chi1Opt = False # Optimize lambda for chi2 = 1
zWeight = 1.0
blockyModel = True

tdip.invertRhoa(mesh=meshInit,lam=lam,robustData=robustData,verbose=verbose,maxIter=1)

ertMgr = tdip.ERT
inv = ertMgr.fw.inv
inv.setModel(pg.Vector(len(tdip.res),np.abs(np.median(tdip.data["rhoa"]))))
inv.setReferenceModel(pg.Vector(len(tdip.res),np.abs(np.median(tdip.data["rhoa"]))))
inv.setRobustData(robustData)
inv.setVerbose(verbose)
inv.setLambda(lam)
inv.setBlockyModel(blockyModel)
inv.setDeltaPhiAbortPercent(1)
inv.setRelativeError(pg.Vector(dataErr))
#inv.setAbsoluteError(dataErr)
tdip.res = inv.run()
tdip.showResistivity(cMap='jet')


data["err"]=deepcopy(np.asarray((dataErr))) # apparently invertRhoa modifies the error column. We reinitialize to the original values here (useful if a trimmed dataset is saved, see after)


# ### Using this inversion, we filter the data and re-run the inversion on a finer mesh (round 2)

# In[70]:


simData=ertMgr.inv.response
obsData=data["rhoa"]
relError = ((abs(obsData-simData))/obsData * 100) # Relative error in % for each data point
cutOffError=5 # Maximum accepted relative error on the reconstruction (%)
plotCompareResults(data, simData)
if inv.relrms() > cutOffError:
    print(f'The inversion run leads to a relative rms of {inv.relrms()} % (chi² = {inv.chi2()}) in {inv.iter()} iterations. \nThis is too high, the dataset will be filtered')
    idxRemove = relError > cutOffError
    idxKeep = [not(i) for i in idxRemove]
    print(f'Out of the initial {len(idxRemove)} data points, {sum(idxKeep)} are kept for inversion.')
    data.remove(idxRemove)
    data.save(directory + 'TrimmedDataset.ohm')
    idxKeep = [not(i) for i in idxRemove]
    relErrorBis = ((abs(data['rhoa']-np.asarray(inv.response())[idxKeep]))/data['rhoa'] * 100)
    plotCompareResults(data, np.asarray(inv.response())[idxKeep])
else:
    print(f'The inversion run leads to a relative rms of {inv.relrms()} % (chi² = {inv.chi2()}) in {inv.iter()} iterations.')
# Re-run the inversion on a finner mesh:

#mgr.setMesh(meshFine)
#pg.show(mgr.paraDomain)

if os.path.exists(directory+"TrimmedDataset.ohm"): # check if TrimmedDataset.ohm exists
    tdip = pb.TDIPdata(directory+"TrimmedDataset.ohm") 
    data = tdip.data
    dataErr=deepcopy(np.asarray(data["err"])) # define the error
else:
    tdip = pb.TDIPdata(directory+dataFile) 
    data=tdip.data
    if not data.allNonZero('k'):
        data['k'] = createGeometricFactors(data, numerical=True)
    if not data.allNonZero('rhoa'):
        data['rhoa'] = data['r']*data['k']
    if not data.allNonZero('err'):
        absoluteError = 0.001
        relativeError = 0.03
        data['err'] = relativeError + absoluteError / data['rhoa'] # In % of the resistance
        dataErr=deepcopy(np.asarray(data["err"])) # define the error
    else:
        print(data["err"])
        dataErr=deepcopy(np.asarray(data["err"])) # define the error
    
mgr = pb.ERTManager(data)      
meshFine = mgr.createMesh(data, paraDepth=paraDepth, quality=quality, paraMaxCellSize=paraMaxCellSize, paraDX=paraDX)
mgr.setMesh(meshFine)

tdip.invertRhoa(mesh=meshFine,lam=lam,robustData=robustData,verbose=verbose,maxIter=1) 
ertMgr = tdip.ERT
inv = ertMgr.fw.inv
inv.setRobustData(robustData)
inv.setVerbose(verbose)
inv.setModel(pg.Vector(len(tdip.res),np.abs(np.median(tdip.data["rhoa"]))))
inv.setReferenceModel(pg.Vector(len(tdip.res),np.abs(np.median(tdip.data["rhoa"]))))
inv.setLambda(lam)
inv.setRelativeError(pg.Vector(dataErr))
inv.setBlockyModel(blockyModel)
inv.setDeltaPhiAbortPercent(1)
tdip.res = inv.run()


# ### Display final resistivity model and fit

# In[71]:

RMSERT=f'The inversion run of ERT data leads to a relative rms of {inv.relrms()} % (chi² = {inv.chi2()}) in {inv.iter()} iterations.'
print(RMSERT)
tdip.coverage = ertMgr.coverage()
tdip.pd = ertMgr.paraDomain
tdip.showResistivity(cMap='jet',cMin=10,cMax=5000,xlabel=xlabel, ylabel=ylabel)
plt.savefig(directory+saveDir+'Resistivity_model.png',transparent=True)
simData=ertMgr.inv.response
obsData=data["rhoa"]
plotCompareResults(data, simData)
plt.savefig(directory+saveDir+'Resistivity_fit.png',transparent=True)


# ### Inversion chargeability with controlled parameters

# In[73]:


ma = tdip.data['ip']*0.001  # chargeability expressed in V/V
    
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
if not data.allNonZero('ipErr'):
    INV.setAbsoluteError(np.ones_like(res)*0.002)
    print('No IP error available, default value of 0.002 V/V used')
else:
    print('IP error available')
    INV.setAbsoluteError(tdip.data['ipErr']*0.001)
        
INV.setLambda(lam)
INV.setRobustData(robustData)
INV.setBlockyModel(blockyModel)
INV.setDeltaPhiAbortPercent(1)

tdip.m = INV.run()
tdip.mafwd = INV.response()

tdip.showChargeability(cMap='jet',cMin = 0, cMax = 25,xlabel=xlabel, ylabel=ylabel)
plt.savefig(directory+saveDir+'Chargeability_model.png',transparent=True)

RMSIP=f'The inversion run of IP data leads to a relative rms of {INV.relrms()} % (chi² = {INV.chi2()}) in {INV.iter()} iterations.'
print(RMSIP)
RMS=[RMSERT,RMSIP]

with open(directory+saveDir+'RMSResults.txt',"w") as f:
    for s in RMS:
        f.write(str(s)+"\n")

# ### Export of results in vtk

# In[74]:


tdip.pd.addData("Resistivity/Ohmm", tdip.res)
tdip.pd.addData("coverage(log10)", tdip.coverage)
tdip.pd.addData("Phase/mrad", tdip.m*1000)
tdip.pd.exportVTK(directory+saveDir+'dcinv.result')

