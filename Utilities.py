import numpy as np

def SaveArrayTXT(filename='DefaultArrayName.txt', array=None):
    if array is not None:
        np.savetxt(filename,array,fmt='%d',delimiter='\t')
