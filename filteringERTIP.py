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
from matplotlib import pyplot as plt

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
        popt, pcov = curve_fit(chargeaOpti, times, IP, np.asarray([10, 0.1, 0]), maxfev=10000, method='trf')
        alpha   = popt[0]
        beta    = popt[1]
        eps     = popt[2]
        if graphs and (value % 100 == 0):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(times, IP, 'or')
            ax.plot(times, chargeaOpti(times, alpha, beta, eps))
            ax.set_title(f'R = {res}')
            plt.show(block=False)
        value += 1
    except:
        nbFails +=1
        alpha   = np.NaN
        beta    = np.NaN
        eps     = np.NaN
        if graphs:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(times, IP, 'or')
            ax.set_title(f'R = {res}')
            plt.show(block=False)
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

## Main parametres:
fileName = './AgatheData/PP3_BerleurCrossFaille.txt'
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
timesIP = (np.cumsum(deltaT) + delayT) * 1000   # Time for the windows (msec)
print(f'Delay Time (sec): {delayT}')
print(f'Window times ({len(deltaT)}): {deltaT}')
print(f'Number of data points: {dataSet.shape[0]} ({dataSet.shape[0]/nbDataInit*100} %)')

ipData = [col for col in dataSet if col.startswith('IP #')]
print(dataSet[ipData].describe())

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


nbBins = 10

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
    stdMisfit[i] = np.std(misfits)
ax.plot(binsCenters, stdMisfit, 'dk')

def errorIP(R, a, b):
    return a * np.power(R, b)

popt, _ = curve_fit(errorIP, binsCenters, stdMisfit, bounds=([0, -np.inf], [np.inf, 0]))
a = popt[0]
b = popt[1]
print(f'Error model for IP: s = {a}R^{b}')
ax.plot(binsCenters, errorIP(binsCenters, a, b), 'k')

def errorIP_r(R, c, d):
    return np.divide(c, R) + d

popt, _ = curve_fit(errorIP_r, binsCenters, stdMisfit, bounds=(0, np.inf))
c = popt[0]
d = popt[1]
print(f'Error model for IP(r): s = {c}/R + {d}')
ax.plot(binsCenters, errorIP_r(binsCenters, c, d), 'g')

plt.show()

## Error model for R:


