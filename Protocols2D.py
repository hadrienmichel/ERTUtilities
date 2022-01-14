import numpy as np
from Utilities import SaveArrayTXT, ArrayToXML, SampleReciprocals

def CreateDDN(nElec:int=64, nMax:int=6, parralelizeMin:int=None):
    '''This function creates a dipole-dipole array for an in-line measurement
    The input parameters are :
        - nElec: the number of electrodes in the line (default = 64)
        - nMax: the maximum n factor (default = 6)
        - parralelizeMin: the minimum number of measures per injection to keep (default=None)

    It outputs a numpy array that contains the array and can be saved in a simple txt file.
    '''
    array = []
    for A in np.arange(start = 1, stop = nElec+1):
        for a in np.arange(start = 1, stop = int(nElec/2)):
            B = A + a
            for n in np.arange(start = 1, stop = nMax+1):
                M = B + n*a 
                N = M + a
                if np.max([A, B, M, N]) <= nElec:
                    array.append([B, A, M, N])
    array = np.asarray(array)
    if parralelizeMin is not None:
        if parralelizeMin > nMax:
            print('Impossible to measure more at least {} dipoles per injection with a maximum n factor of {}.'.format(parralelizeMin, nMax))
            parralelizeMin = nMax
        uniqueAB, countsAB = np.unique(array[:,:2], axis=0, return_counts=True)
        valKeep = uniqueAB[countsAB==parralelizeMin]
        array = []
        for AB in valKeep:
            A, B = AB[0], AB[1]
            a = B - A
            for n in np.arange(start=1, stop = parralelizeMin+1):
                M = B + n*a 
                N = M + a 
                array.append([B, A, M, N])
        array = np.asarray(array)
    return array

def CreateGradient(nElec:int=64, s:int=7):
    '''This function creates a gradient array for an in-line measurement
    The input parameters are :
        - nElec: the number of electrodes in the line (default = 64)
        - s: the s factor (default = 7)
        - parralelizeMin: the minimum number of measures per injection to keep (default=None)

    It outputs a numpy array that contains the array and can be saved in a simple txt file.
    '''
    array = []
    for A in np.arange(start = 1, stop = nElec+1):
        for a in np.arange(start = 1, stop = int(nElec/2)):
            B = A + (s+2)*a 
            for M in np.arange(start = A+a, stop = A + (s+1)*a, step = a):
                N = M + a
                if np.max([A, B, M, N]) <= nElec:
                    array.append([A, B, M, N])
    array = np.asarray(array)
    return array

if __name__=='__main__':
    n = 6
    nElec = 64
    DDNArray = CreateDDN(nMax = n, nElec = nElec, parralelizeMin=6)
    print('Default DDN6 array (length = {}):'.format(len(DDNArray)))
    #print(DDNArray)
    ArrayToXML(DDNArray)
    DDNArrayReciprocals = SampleReciprocals(DDNArray)
    print('The corresponding DDN6 sampled reciprocal array (length = {}):'.format(len(DDNArrayReciprocals)))
    ArrayToXML(DDNArrayReciprocals)
    s = 7
    GradientArray = CreateGradient(s = s, nElec = nElec)
    print('Default Gradient (s=7) array (length = {})'.format(len(GradientArray)))
    #print(GradientArray)
    ArrayToXML(GradientArray)
    GradientArrayReciprocals = SampleReciprocals(GradientArray)
    print('The corresponding Gradient (s=7) sampled reciprocal array (length = {}):'.format(len(GradientArrayReciprocals)))
    ArrayToXML(GradientArrayReciprocals)