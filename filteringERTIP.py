'''Filtering and modelisation of IP and ERT data and error estimation.

Adapted from:
Flores Orozco A., Gallistl J., Bücker M. & Williams K. (2018). Decay curve analysis for
data error quantification in time-domain induced polarization imaging. GEOPHYSICS, 83 (2),
1MA-Z8. https://doi.org/10.1190/geo2016-0714.1

Use only for IP filtering or when reciprocals are not present.
'''
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize_scalar
from math import isnan
from matplotlib import pyplot as plt
from pygimli.frameworks import fit as fitFunc

import ExtractLines as el

import re
# From StackOverflow: https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string
numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
floatsExtractor = re.compile(numeric_const_pattern, re.VERBOSE)

# Counters initializer:
nbFails = 0
value = 0
# Global values:
graphs = False

electrodesMatrix = []

## Defining common functions:
def computeChargeability(IP, windows):
    a = 0
    for i, ti in enumerate(windows):
        a += IP[i]*ti
    return a/np.sum(windows)

def chargeaOpti(t, alpha, beta, eps):
    return alpha * np.power(t, -beta) + eps

def fitChargea(IP, times, res):
    global nbFails, value, graphs
    try:
        # popt, _ = fitFunc(chargeaOpti, data=IP, t=times, startModel=[10.0, 0.1, 0.1], verbose=True)
        popt, _ = curve_fit(chargeaOpti, times, IP, np.asarray([10, 0.1, 0]), method='trf', maxfev=2000, loss='soft_l1')
        alpha   = popt[0]
        beta    = popt[1]
        eps     = popt[2]
        value += 1
        if (graphs and (value % 100 == 0)) or (graphs and res < 0.5):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(times, IP, 'or')
            ax.plot(times, chargeaOpti(times, alpha, beta, eps))
            ax.set_title(f'curve nb. {value} - R = {res}')
            plt.show(block=True)
    except:
        value += 1
        nbFails +=1
        alpha   = np.NaN
        beta    = np.NaN
        eps     = np.NaN
        if graphs:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(times, IP, 'or')
            ax.set_title(f'curve nb. {value} - R = {res}')
            plt.show(block=True)
    return alpha, beta, eps

def computeRMSD(IP, times, alpha, beta, eps):
    return np.sqrt(1/len(IP) * np.sum(np.power(chargeaOpti(times, alpha, beta, eps) - IP, 2)))

def referenceMedian(row, df, times):
    # Find the index of similar A and B:
    currDf = df[(df['A(x)'] == row['A(x)']) & (df['A(y)'] == row['A(y)']) & (df['A(z)'] == row['A(z)']) &
                (df['B(x)'] == row['B(x)']) & (df['B(y)'] == row['B(y)']) & (df['B(z)'] == row['B(z)'])] # Select the rows with the same a and b as the idxRow one.
    # Compute the reference:
    ref = np.median(np.asarray([chargeaOpti(times, rowLoop['alpha'], rowLoop['beta'], rowLoop['eps']) for _, rowLoop in currDf.iterrows()]), axis=0)
    return ref

def rmsdRef(ref, fit):
    return np.sqrt(1/len(ref) * np.sum(np.power(ref - fit, 2)))

def computeShift(ref, fit):
    ku = np.nan
    kd = np.nan
    def misfitShift(k):
        return rmsdRef(ref, fit+k)
    kOpt = minimize_scalar(misfitShift)
    if kOpt.success:
        kOpt = kOpt.x
        if kOpt < 0:
            kd = -kOpt
        else:
            ku = kOpt
    else:
        print('Failed to minimize the functional')
    return ku, kd

def assignElecId(row, df):
    global electrodesMatrix
    keyCol = ['aId', 'bId', 'mId', 'nId']
    corrKey = [['A(x)', 'A(y)', 'A(z)'],
               ['B(x)', 'B(y)', 'B(z)'],
               ['M(x)', 'M(y)', 'M(z)'],
               ['N(x)', 'N(y)', 'N(z)']]
    for i, key in enumerate(keyCol):
        if isnan(row[key]):
            ids = df[(df[corrKey[i][0]] == row[corrKey[i][0]]) & 
                     (df[corrKey[i][1]] == row[corrKey[i][1]]) &
                     (df[corrKey[i][2]] == row[corrKey[i][2]])].index.tolist()
            try:
                electrodeId = electrodesMatrix.index(row[corrKey[i]].tolist())
            except:
                # Electrode not in array
                electrodesMatrix.append(row[corrKey[i]].tolist())
                electrodeId = electrodesMatrix.index(row[corrKey[i]].tolist())
            df[key][ids] = electrodeId



## Main parametres:
fileName = './AgatheData/PP9_ChessionAmont.txt'
errorModelling = True

## Reading the file using pandas
deltaT = None
f = open(fileName, 'r')
for line in f.readlines():
    if deltaT is not None:
        break
    if 'IP_WindowSecList' in line:
        deltaT = [float(i) for i in floatsExtractor.findall(line)]
lineHeaders, lineEnd = el.getLines(filename=fileName)
dataSet = pd.read_csv(fileName, skiprows=lineHeaders, nrows=lineEnd-lineHeaders-2, sep='\t', index_col=False) # Need to check for 2 rows too much. Is it the case for every files???
nbDataInit = dataSet.shape[0]
delayT = deltaT[0]
deltaT = deltaT[1:]
print(f'Delay Time (sec): {delayT}')
print(f'Window times ({len(deltaT)}): {deltaT}')
print(f'Number of data points: {dataSet.shape[0]} ({dataSet.shape[0]/nbDataInit*100} %)')

elecA = np.unique(dataSet[['A(x)', 'A(y)', 'A(z)']].to_numpy(), axis=0)
elecB = np.unique(dataSet[['B(x)', 'B(y)', 'B(z)']].to_numpy(), axis=0)
elecM = np.unique(dataSet[['M(x)', 'M(y)', 'M(z)']].to_numpy(), axis=0)
elecN = np.unique(dataSet[['N(x)', 'N(y)', 'N(z)']].to_numpy(), axis=0)
electrodesMatrix = np.unique(np.vstack((elecA, elecB, elecM, elecN)), axis=0).tolist()

ipData = [col for col in dataSet if col.startswith('IP #')]
print(dataSet[ipData].describe())
while (dataSet[ipData[0]] == 0).all():
    # The first column of the ip dataset was not registered
    # Merging with the second column
    print('The first window was not aquiered.\nMerging with the second column!')
    # Swithc all data one colum:
    dataSet[ipData[:-1]] = dataSet[ipData[1:]]
    # Remove the last colum:
    deltaT[1] = deltaT[0] + deltaT[1]
    deltaT = deltaT[1:]
    dataSet.drop(columns=ipData[-1])
    ipData = ipData[:-1]
print(dataSet[ipData].describe())
tmp = np.cumsum(np.hstack([[delayT], deltaT]))

timesIP = np.divide(tmp[:-1]+tmp[1:], 2) # (np.cumsum(deltaT) + delayT)#  * 1000   # Time for the windows (msec)
print(timesIP)
## Fitting of a power law mf(t) = alpha*t^(-beta) + eta:
dataSet['Chargeability (mV/V)'] = dataSet.apply(lambda row: computeChargeability(row[ipData], deltaT), axis=1)
dataSet[['alpha', 'beta', 'eps']] = dataSet.apply(lambda row: fitChargea(row[ipData], timesIP, row['R(Ohm)']), axis=1, result_type='expand')
dataSet = dataSet.dropna(subset=['alpha', 'beta', 'eps'])
print(dataSet[['Chargeability (mV/V)', 'alpha', 'beta']].describe())
print(f'Number of failed fit: {nbFails}')
if graphs: 
    plt.show()

## Fisrt filter (non decaying curves) - alpha and beta of diferent signs:
prod = np.multiply(dataSet['alpha'], dataSet['beta'])
dataSet = dataSet[prod>0]
print(f'Number of data points after the first filter: {dataSet.shape[0]} ({dataSet.shape[0]/nbDataInit*100} %)')

## Root mean square deviation computation:
dataSet['rmsd'] = dataSet.apply(lambda row: computeRMSD(row[ipData], timesIP, row['alpha'], row['beta'], row['eps']), axis=1)

## Construction of reference curves for each AB couples:
dataSet['Reference'] = dataSet.apply(lambda row: referenceMedian(row, dataSet, timesIP), axis=1) # Median curve
dataSet['Fitted'] = dataSet.apply(lambda row: chargeaOpti(timesIP, row['alpha'], row['beta'], row['eps']), axis=1) # Fitted curve
dataSet['Measured'] = dataSet.apply(lambda row: np.asarray(row[ipData]).flatten(), axis=1) # Measured curve

dataSet[['ku', 'kd']] = dataSet.apply(lambda row: computeShift(row['Reference'], row['Fitted']), axis=1, result_type='expand')

## Second filter (far from reference curves):
stdKu = dataSet['ku'].std()
stdKd = dataSet['kd'].std()
# Calculating the noise level:
kuMult, kdMult = 3, 3 # General case
if 3*stdKu >  2*dataSet['Chargeability (mV/V)'].median():
    print('The dataset is noisy!')
    # Noisy dataset
    kuMult, kdMult = 1.5, 1
elif 3*stdKu < dataSet['Chargeability (mV/V)'].median():
    # Clean dataset
    print('The dataset is clean!')
    kuMult, kdMult = 4, 4

dataSet = dataSet[(dataSet['ku'] < kuMult*stdKu) | (dataSet['kd'] < kdMult*stdKd)]
print(f'Number of data points after the second filter: {dataSet.shape[0]} ({dataSet.shape[0]/nbDataInit*100} %)')
print(dataSet['Chargeability (mV/V)'].describe())
dataSet = dataSet.reset_index(drop=True)
## Standard deviation estimation for IP:
dataSet['DCAmisfit'] = dataSet.apply(lambda row: row['Fitted']-row['Measured'], axis=1)


nbBins = 5

fig = plt.figure()
ax = fig.add_subplot(111)
for _, row in dataSet.iterrows():
    ax.plot(np.ones((len(timesIP),))*np.abs(row['R(Ohm)']), row['DCAmisfit'], 'or')
ax.set_xscale('log')
ax.set_xlabel('R (Ohm)')
ax.set_ylabel('Misfit (in mV/V)')
ax.set_ylim([-25, 25])
ax.grid()

bins = np.logspace(np.log10(dataSet['R(Ohm)'].min()), np.log10(dataSet['R(Ohm)'].max()), num=nbBins+1)
binsCenters = (bins[1:] + bins[:-1])/2
stdMisfit = np.zeros((nbBins,))
for i in range(nbBins):
    currData = dataSet[(dataSet['R(Ohm)']>=bins[i]) & (dataSet['R(Ohm)']<bins[i+1])]
    misfits = np.array(currData['DCAmisfit'].to_list()).flatten()
    stdMisfit[i] = np.std(misfits) # np.mean(np.absolute(misfits - np.mean(misfits))) #  
ax.plot(binsCenters, stdMisfit, 'dk')

def errorIP(R, a, b):
    return a * np.power(R, b)

popt, _ = curve_fit(errorIP, binsCenters, stdMisfit, bounds=([0, -10], [10, 0]), method='trf', loss='soft_l1')
a = popt[0]
b = popt[1]
print(f'Error model for IP: s = {a}R^{b}')
ax.plot(binsCenters, errorIP(binsCenters, a, b), 'k')

dataSet['ipErr'] = dataSet.apply(lambda row: errorIP(np.abs(row['R(Ohm)']), a, b), axis=1)

def errorIP_r(R, c, d):
    return np.divide(c, R) + d

popt, _ = curve_fit(errorIP_r, binsCenters, stdMisfit, bounds=([0, 0], [10, 10]), method='trf', loss='soft_l1')
c = popt[0]
d = popt[1]
print(f'Error model for IP(r): s = {c}/R + {d}')
ax.plot(binsCenters, errorIP_r(binsCenters, c, d), 'g')

## Error model for R:

dataSet['err'] = dataSet.apply(lambda row: (c + d*row['R(Ohm)'])/row['R(Ohm)'], axis=1)

dataSet.hist(column=['Var(%)', 'err', 'ipErr', 'Chargeability (mV/V)'])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataSet['Chargeability (mV/V)'], dataSet['ipErr'],'.b')
ax.set_ylabel('Error (mV/V)')
ax.set_xlabel('Chargeability (mV/V)')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataSet['R(Ohm)'], dataSet['err'],'.b', label='Reciprocals')
ax.plot(dataSet['R(Ohm)'], dataSet['Var(%)'], '.r', label='Repetition')
ax.set_ylabel('Error (%)')
ax.set_xlabel('R (Ohm)')
plt.legend()

plt.show()

## Saving dataset in BERT formatting:

# 1) Ask for the file to save
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

file_path = filedialog.asksaveasfilename(defaultextension='.ohm', filetypes=(('Resistivity file', '.ohm'), ('Data file', '.dat'), ('Unified data format', '.udf')), title='Save data as . . .')

# 2) Sort the data:
dataSet = dataSet.assign(aId=np.NaN, bId=np.NaN, mId=np.NaN, nId=np.NaN)
dataSet['elecId'] = dataSet.apply(lambda row: assignElecId(row, dataSet), axis=1) #pd.read_csv(electrodesFile, sep='\t', index_col=False)
print(f'The number of electrodes in the dataset is: {len(electrodesMatrix)}')
# elecCorr = np.asarray(dfElec[['Xequi', 'Yequi', 'Zequi']])
# if 'Xchange' in dfElec.columns: 
#     # The electrodes position is changed (added topography for example)
#     elecChanges = np.asarray(dfElec[['Xchange', 'Ychange', 'Zchange']])
# else:
#     elecChanges = elecCorr
# for el in range(len(electrodes)):
#     idx = np.where((elecCorr == electrodes[el,:-1]).all(axis=1))
#     electrodes[el,:-1] = elecChanges[idx,:]
# dataset = np.asarray(dataSet[['A(adr)','B(adr)','M(adr)','N(adr)','R(Ohm)', 'Var(%)']])
# tol = 1e-3 #Tolerance on the elctrode position ~ cm
# for data in range(len(dataset)):
#     aIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,0])<tol)[0][0] + 1
#     bIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,1])<tol)[0][0] + 1
#     mIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,2])<tol)[0][0] + 1
#     nIdx = np.where(np.abs(electrodes[:,-1]-dataset[data,3])<tol)[0][0] + 1
#     dataset[data,:4] = [aIdx, bIdx, mIdx, nIdx]
dataset = np.asarray(dataSet[['aId','bId','mId','nId','R(Ohm)','err','Chargeability (mV/V)','ipErr']])
datasetElectrodes = (dataset[:,:4] + 1).astype(np.int)
dataset = dataset[:,4:]
electrodes = np.asarray(electrodesMatrix)
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
    file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(*datasetElectrodes[i,:], *dataset[i,:]))
file.close()
print('File saved succesfully!')
