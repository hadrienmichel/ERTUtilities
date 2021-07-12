import numpy as np

def SaveArrayTXT(filename='DefaultArrayName.txt', array=None):
    if array is not None:
        np.savetxt(filename,array,fmt='%d',delimiter='\t')
    
def SortArray(array=None, type='ABMN'):
    if array is not None:
        if type == 'ABMN':
            array = array[array[:,3].argsort()]
            array = array[array[:,2].argsort(kind='mergesort')]
            array = array[array[:,1].argsort(kind='mergesort')]
            array = array[array[:,0].argsort(kind='mergesort')]
        if type == 'NMBA':
            array = array[array[:,0].argsort()]
            array = array[array[:,1].argsort(kind='mergesort')]
            array = array[array[:,2].argsort(kind='mergesort')]
            array = array[array[:,3].argsort(kind='mergesort')]
    return array

def CreateReciprocals(array):
    '''Create the reciprocals array for a given array in input.
    '''
    if array is not None:
        return array[:,[2, 3, 0, 1]]

def SampleReciprocals(array, sampling=0.1):
    '''Sample the reciprocals array for the injections
    '''
    if array is not None:
        reciprocals = CreateReciprocals(array)
        ABInj = np.unique(reciprocals[:,:2], axis=0)
        nbSamples = int(len(ABInj)*sampling)
        print(f'{nbSamples} injection in reciprocals')
        idxKeep = np.random.choice(np.arange(len(ABInj)),nbSamples,replace=False)
        ABInjKeep = ABInj[idxKeep,:]
        ReciArray = np.asarray([0, 0, 0, 0])
        for AB in ABInjKeep:
            # Find the corresponding arrays:
            A, B = AB[0], AB[1]
            idxCurr = np.where((reciprocals[:,:2] == (A, B)).all(axis=1))
            ReciArray = np.vstack([ReciArray, np.squeeze(reciprocals[idxCurr,:])])
        ReciArray = ReciArray[1:,:]
        return ReciArray

