import numpy as np
from numpy.core import multiarray
from Protocols2D import CreateGradient, CreateDDN
from Utilities import SaveArrayTXT
from copy import deepcopy

def CreateMultiLineGradient(nLines:int=4, nElecLine:int=32, s:int=7):
    simpleGradient = CreateGradient(nElec = nElecLine, s = s)
    multiGradient = np.asarray([0, 0, 0, 0])
    for lineInj in range(nLines):
        for lineMeas in np.arange(start=lineInj, stop=nLines):
            array = deepcopy(simpleGradient)
            array[:,:2] += lineInj*nElecLine
            array[:,2:] += lineMeas*nElecLine
            multiGradient = np.vstack([multiGradient,array])
    multiGradient = multiGradient[1:,:]
    return multiGradient

def CreateMultiLineDDN(nLines:int=4, nElecLine:int=32, nMax:int=6):
    simpleDDN = CreateDDN(nElec = nElecLine, nMax = nMax, parralelizeMin = nMax)
    multiDDN = np.asarray([0, 0, 0, 0])
    for lineInj in range(nLines):
        for lineMeas in np.arange(start=lineInj, stop=nLines):
            array = deepcopy(simpleDDN)
            array[:,:2] += lineInj*nElecLine
            array[:,2:] += lineMeas*nElecLine
            multiDDN = np.vstack([multiDDN,array])
    multiDDN = multiDDN[1:,:]
    return multiDDN

if __name__ == "__main__":
    multiDDN = CreateMultiLineGradient()
    print('Default 4 lines Gradient array (length = {}):'.format(len(multiDDN)))
    print(multiDDN)
    multiDDN = CreateMultiLineDDN()
    print('Default 4 lines DDN6 array (length = {}):'.format(len(multiDDN)))
    print(multiDDN)

