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
