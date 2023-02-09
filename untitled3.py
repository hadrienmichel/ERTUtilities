# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 14:52:00 2021

@author: tomde
"""
import numpy as np
from Utilities import SaveArrayTXT 

from Protocols2D import CreateDDN

def CreateGradient(nElec:int=32, s:int=7):
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

output=CreateGradient()


SaveArrayTXT(filename='GD7.txt', array=output)
print(output)