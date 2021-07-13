import numpy as np
from Protocols2D import CreateGradient, CreateDDN
from Utilities import ArrayToXML, SampleReciprocals, SaveArrayTXT, SortArray
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
    return SortArray(multiGradient)

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
    return SortArray(multiDDN)

def CreateTom(testCase:int=0):
    nElecLine = 32
    s = 7
    simpleGradient = CreateGradient(nElec=nElecLine, s=s)
    ABInj = np.unique(simpleGradient[:,:2],axis=0)
    if testCase == 0:
        # Case with full 3D array:
        return CreateMultiLineGradient()
    elif testCase == 1:
        # Case with 2 injection lines in the center
        # We only inject on the middle lines and measure fully on those.
        # Cross line measurements are performed sparsly on the other lines
        FullArray = np.asarray([0, 0, 0, 0])
        injectionLines = [1, 2]
        crossLines = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
        for lineInj in injectionLines:
            array = deepcopy(simpleGradient)
            array[:,:2] += lineInj*nElecLine
            array[:,2:] += lineInj*nElecLine
            FullArray = np.vstack([FullArray, array])
            for AB in ABInj:
                A = AB[0]
                B = AB[1]
                a = int((B-A)/(s+2))
                M = A + 4*a #int((B-A)/2 - a/2)
                N = M + a #int((B-A)/2 + a/2)
                for lineMeas in crossLines[lineInj]:
                    ACurr = A + lineInj*nElecLine
                    BCurr = B + lineInj*nElecLine
                    MCurr = M + lineMeas*nElecLine
                    NCurr = N + lineMeas*nElecLine
                    FullArray = np.vstack([FullArray, [ACurr, BCurr, MCurr, NCurr]])
        FullArray = FullArray[1:,:]
    elif testCase==2:
        # Case with injection on a diagonal on outer lines and measures cross-lined in the middle
        simpleGradient10 = CreateGradient(nElec= nElecLine, s=10)
        # From left to right:
        array1 = deepcopy(simpleGradient10)
        array1[:,1] += 3*nElecLine
        array1[:,2] += 1*nElecLine
        array1[:,3] += 2*nElecLine
        # From right to left:
        array2 = deepcopy(simpleGradient10)
        array2[:,0] += 3*nElecLine
        array2[:,2] += 2*nElecLine
        array2[:,3] += 1*nElecLine

        FullArray = np.vstack([array1, array2])

    return SortArray(FullArray)

if __name__ == "__main__":
    multiGradient = CreateMultiLineGradient()
    print('Default 4 lines Gradient array (length = {}):'.format(len(multiGradient)))
    # print(multiGradient)
    multiDDN = CreateMultiLineDDN()
    print('Default 4 lines DDN6 array (length = {}):'.format(len(multiDDN)))
    # print(multiDDN)
    # Testing Tom's protocols:
    multiTom1 = CreateTom(testCase=1)
    print('Multiple methods array (Tom, case=1) (length = {} - {} injections):'.format(len(multiTom1), len(np.unique(multiTom1[:,:2], axis=0))))
    # print(multiTom1)
    multiTom2 = CreateTom(testCase=2)
    print('Multiple methods array (Tom, case=2) (length = {} - {} injections):'.format(len(multiTom2), len(np.unique(multiTom2[:,:2], axis=0))))
    # print(multiTom2)
    multiFullTom = np.vstack([multiTom1, multiTom2])
    print('Multiple methods array (Tom, case=3) (length = {}):'.format(len(multiFullTom)))
    # print(multiFullTom)
    print('{} injections for the full array!'.format(len(np.unique(multiFullTom[:,:2], axis=0))))
    # SaveArrayTXT('MultiTomFull.txt',multiFullTom)
    RecipMultiFullTom = SampleReciprocals(multiFullTom, sampling=0.15)
    print('Multiple methods array (Tom, case=3) (length = {}):'.format(len(RecipMultiFullTom)))
    # print(RecipMultiFullTom)
    # SaveArrayTXT('MultiTomFullReciprocals.txt',RecipMultiFullTom)
    ArrayToXML(multiFullTom)
    ArrayToXML(RecipMultiFullTom)
    print('End of File . . .')
