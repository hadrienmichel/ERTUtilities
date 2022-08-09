import numpy as np
from Utilities import SaveArrayTXT, ArrayToXML, SampleReciprocals, SortArray

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
                    array.append([A, B, N, M])
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
                array.append([A, B, N, M])
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

def makeBorehole(array, nElecSurface:int=64, posBorehole:int=32.5, nElecBorehole:int=32):
    addArray = []
    for abmn in array:
        A, B, M, N = abmn[0], abmn[1], abmn[2], abmn[3]
        if (M > posBorehole) or (N > posBorehole):
            Madd, Nadd = M, N
            if M > posBorehole: Madd = M + np.floor(posBorehole)
            if N > posBorehole: Nadd = N + np.floor(posBorehole)
            addArray.append([A, B, Madd, Nadd])
            if (A+B)/2 != posBorehole:
                Aadd, Badd = A, B
                if M < posBorehole: Madd = nElecSurface+1 - M
                if N < posBorehole: Nadd = nElecSurface+1 - N
                if A < posBorehole: Aadd = nElecSurface+1 - A
                if B < posBorehole: Badd = nElecSurface+1 - B
                addArray.append([Aadd, Badd, Madd, Nadd])
    addArray = np.asarray(addArray)
    array = SortArray(np.vstack((array, addArray)).astype(int))
    return array


if __name__=='__main__':
    n = 6
    nElec = 64
    DDNArray = CreateDDN(nMax = n, nElec = nElec, parralelizeMin=6)
    print('Default DDN6 array (length = {}):'.format(len(DDNArray)))
    print(DDNArray)
    DDN6Bore = makeBorehole(DDNArray)
    print('Default DDN6 with borehole array (length = {}):'.format(len(DDN6Bore)))
    print(DDN6Bore)
    SaveArrayTXT('DDN6Borehole.txt', DDN6Bore)
    DDN6BoreRecip = SampleReciprocals(DDN6Bore)
    SaveArrayTXT('DDN6BoreholeReciprocals.txt', DDN6BoreRecip)

    s = 7
    GradientArray = CreateGradient(s = s, nElec = nElec)
    print('Default Gradient (s=7) array (length = {})'.format(len(GradientArray)))
    GradientBore = makeBorehole(GradientArray)
    print('Default DDN6 with borehole array (length = {}):'.format(len(GradientBore)))
    print(GradientBore)
    SaveArrayTXT('Gradient7Borehole.txt', GradientBore)
    GradientBoreRecip = SampleReciprocals(GradientBore)
    SaveArrayTXT('GradientBoreholeReciprocals.txt', GradientBoreRecip)
    # ArrayToXML(DDNArray)
    # DDNArrayReciprocals = SampleReciprocals(DDNArray)
    # print('The corresponding DDN6 sampled reciprocal array (length = {}):'.format(len(DDNArrayReciprocals)))
    # ArrayToXML(DDNArrayReciprocals)
    # s = 7
    # GradientArray = CreateGradient(s = s, nElec = nElec)
    # print('Default Gradient (s=7) array (length = {})'.format(len(GradientArray)))
    # #print(GradientArray)
    # ArrayToXML(GradientArray)
    # GradientArrayReciprocals = SampleReciprocals(GradientArray)
    # print('The corresponding Gradient (s=7) sampled reciprocal array (length = {}):'.format(len(GradientArrayReciprocals)))
    # ArrayToXML(GradientArrayReciprocals)