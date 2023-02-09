import ExtractLines as el
import numpy as np
import pandas as pd
import math as mt
from matplotlib import pyplot as plt

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Retreived and modified from : https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return np.greater(modified_z_score,thresh)

normalMeas = './VSL_DATA/Project2_G7_ABMN_3.txt'#'./VSL_Data/Project2_DDN6_1.txt' #'./VSL_DATA/Project2_G7_ABMN_3.txt'
recipMeas = './VSL_Data/Project2_G7_RECIP_2.txt'#'./VSL_Data/Project2_DDN6_R_S_1.txt'#'./VSL_Data/Project2_G7_RECIP_2.txt'
electrodesFile = './VSL_Data/Line2_Electrodes.txt'
varMax = 0.5

## Load the datasets:
# 1) Normal measurements
lineHeaders, lineEnd = el.getLines(filename=normalMeas)
normal = pd.read_csv(normalMeas, skiprows=lineHeaders, nrows=lineEnd-lineHeaders, sep='\t', index_col=False)
normal = normal[normal['Var(%)']<varMax]
print(normal)
# Reciprocal measurements
lineHeaders, lineEnd = el.getLines(filename=recipMeas)
reciprocal = pd.read_csv(recipMeas, skiprows=lineHeaders, nrows=lineEnd-lineHeaders, sep='\t', index_col=False)
reciprocal = reciprocal[reciprocal['Var(%)']<varMax]
print(reciprocal)

# # To remove for other datasets (change A and B):
# normal[['A(x)','A(y)','A(z)','A(adr)']], normal[['B(x)','B(y)','B(z)','B(adr)']] = normal[['B(x)','B(y)','B(z)','B(adr)']], normal[['A(x)','A(y)','A(z)','A(adr)']]
# reciprocal[['M(x)','M(y)','M(z)','M(adr)']], reciprocal[['N(x)','N(y)','N(z)','N(adr)']] = reciprocal[['N(x)','N(y)','N(z)','N(adr)']], reciprocal[['M(x)','M(y)','M(z)','M(adr)']]
# normal['R(Ohm)'] = -normal['R(Ohm)']
# reciprocal['R(Ohm)'] = -reciprocal['R(Ohm)']
# # End substitution A and B

# Get a list of the electrodes positions:
positionsA = normal[['A(x)','A(y)','A(z)','A(adr)']]
positionsB = normal[['B(x)','B(y)','B(z)','B(adr)']]
positionsM = normal[['M(x)','M(y)','M(z)','M(adr)']]
positionsN = normal[['N(x)','N(y)','N(z)','N(adr)']]
positionsTotal = np.concatenate((positionsA, positionsB, positionsM, positionsN), axis=0)
electrodes = np.unique(positionsTotal,axis=0)
print(electrodes)
print('They are {} electrodes in the dataset.'.format(len(electrodes)))

# normal[['R(Ohm)','Rho-a(Ohm-m)']].hist()
# reciprocal[['R(Ohm)','Rho-a(Ohm-m)']].hist()
# plt.show()

## Buidling the slatter model
# 1) Find the measurements that are done in both arrays:
arrayRecip = np.asarray(reciprocal[['M(adr)', 'N(adr)', 'A(adr)', 'B(adr)']])
Rrecip = np.asarray(reciprocal['R(Ohm)'])
arrayNormal = np.asarray(normal[['A(adr)', 'B(adr)', 'M(adr)', 'N(adr)']])
Rnormal = np.asarray(normal['R(Ohm)'])
RSlater = np.zeros((len(Rrecip),))
eSlater = np.zeros((len(Rrecip),))
for idxRecip, meas in enumerate(arrayRecip):
    idxNorm = (arrayNormal == meas).all(axis=1).nonzero()
    RSlater[idxRecip] = (Rnormal[idxNorm] + Rrecip[idxRecip])/2
    eSlater[idxRecip] = np.abs(Rnormal[idxNorm] - Rrecip[idxRecip])

# 2) Remove outliers (optional):
removeOutliers = True
if removeOutliers:
    outliers = is_outlier(eSlater)
    eSlaterOutliers = eSlater[outliers]
    RSlaterOutliers = RSlater[outliers]
    eSlater = eSlater[np.logical_not(outliers)]
    RSlater = RSlater[np.logical_not(outliers)]

# 3) Find the law of error:
# err = a + b*R, passing by the absolute maximum and the y-axis at a given place
eMaxSlater = np.amax(eSlater)
RMaxSlater = RSlater[eSlater==eMaxSlater]

def errorModelling(a, b, R):
    return np.multiply(R, b) + a

a0 = 0.0
a1 = eMaxSlater
for i in range(20):
    a2 = (a0 + a1) / 2
    e2 = errorModelling(a2, (eMaxSlater-a2)/RMaxSlater, RSlater[eSlater!=eMaxSlater]) - eSlater[eSlater!=eMaxSlater] # We do not take into account the max value as we are sure to pass by this
    maxDistSigned = np.amin(e2)
    if maxDistSigned < 0:
        # Value of a is too low
        a0 = a2
    elif maxDistSigned > 0:
        # Value of a is too high
        a1 = a2
    else:
        # Value of a is exact
        break

aSlater = a2
bSlater = (eMaxSlater-a2)/RMaxSlater

# 4) Show graph with the reciprocal error:
plt.figure()
plt.scatter(RSlater,eSlater,color='b')
if removeOutliers:
    plt.scatter(RSlaterOutliers, eSlaterOutliers, color=[0.5, 0.5, 0.5])
rExample = np.linspace(0.0, 10.0, 100)
plt.plot(rExample, errorModelling(aSlater, bSlater, rExample),color='r')
plt.grid(True)
plt.xlim(left=0.0)
plt.ylim(bottom=0.0)
plt.title('Slater Model of reciprocal error')
plt.xlabel('Resistance [Ohm]')
plt.ylabel('Reciprocal error [Ohm]')
plt.show()

## Save to the *.udf format (BERT):
# 1) Ask for the file to save
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

file_path = filedialog.asksaveasfilename(defaultextension='.dat', filetypes=(('Data file', '.dat'), ('Unified data format', '.udf')), title='Save data as . . .')

# 2) Sort the data:
dfElec = pd.read_csv(electrodesFile, sep='\t', index_col=False)
elecCorr = np.asarray(dfElec[['Xequi', 'Yequi', 'Zequi']])
if 'Xchange' in dfElec.columns: 
    # The electrodes position is changed (added topography for example)
    elecChanges = np.asarray(dfElec[['Xchange', 'Ychange', 'Zchange']])
else:
    elecChanges = elecCorr
for el in range(len(electrodes)):
    idx = np.where((elecCorr == electrodes[el,:-1]).all(axis=1))
    electrodes[el,:-1] = elecChanges[idx,:]
dataset = np.asarray(normal[['A(adr)','B(adr)','M(adr)','N(adr)','R(Ohm)', 'Var(%)']])
tol = 1e-3 #Tolerance on the elctrode position ~ cm
for data in range(len(dataset)):
    aIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,0])<tol)[0][0] + 1
    bIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,1])<tol)[0][0] + 1
    mIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,2])<tol)[0][0] + 1
    nIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,3])<tol)[0][0] + 1
    dataset[data,:-2] = [aIdx, bIdx, mIdx, nIdx]
    dataset[data,-1] = (aSlater + bSlater * dataset[data, -2])/dataset[data, -2]

# 3) Write to the file:
if len(np.unique(electrodes[:,1])) > 1:
    print('The dataset is 3D.')
    is2D = False
else:
    print('The dataset is 2D.')
    is2D = True
file = open(file_path, 'w')
file.write('{}# Number of electrodes\n'.format(len(electrodes)))
if is2D:
    file.write('# x z\n')
else:
    file.write('# x y z\n')
for i in range(len(electrodes)):
    if is2D:
        file.write('{}\t{}\n'.format(electrodes[i,0], electrodes[i,2])) # x and z only
    else:
        file.write('{}\t{}\t{}\n'.format(electrodes[i,0], electrodes[i,1], electrodes[i,2])) # x, y and z
file.write('{}# Number of data\n'.format(len(dataset)))
file.write('# a b m n rho err ip ipErr\n')
for i in range(len(dataset)):
    file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(*dataset[i,:]))
file.close()
