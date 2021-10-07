'''This code enables the rapid extraction of single lines measurements out of a 3D regular array.
The lines can be both along the x and y axis.

Please, complete the parameters in the section below to be able to run this code.
'''
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join

### Input parameters:
# Lines to extract:
yLines = [11, 16, 21] # Lines to be extracted (x-direction)
nameLinesY = ['Line1','Line2','Line3']
xLines = [24] # Lines to be extracted (y-direction)
nameLinesX = ['Line4']
# Directory where the files are located (only .txt files WITH DATA, others will trigger an error!)
directory = './TestFiles/'
# Coordinates of electrodes to remove from the dataset.
coordRemove = [[25.5, 16],
               [25.5, 21],
               [24, 0.5],
               [24, 29.5],
               [24, 30.5],
               [24, 31.5]]
# Maximum variability on the repetition error from the ABEM measurements
varMax = 5.0

def extractLine(db, y=16):
    return db[np.asarray([db['A(y)']==y, db['B(y)']==y, db['M(y)']==y, db['N(y)']==y]).all(axis=0)]

def extractTransect(db, x=24):
    return db[np.asarray([db['A(x)']==x, db['B(x)']==x, db['M(x)']==x, db['N(x)']==x]).all(axis=0)]

def removeFaultyElectrode(db, coordRemove):
    '''Removes the faulty electrodes from the dataset. 
    The positions of the electrodes are stored in coordRemove.
    '''
    nbPoints = len(db)
    for elec in coordRemove:
        db = db[np.asarray([db['A(x)']!=elec[0], db['A(y)']!=elec[1]]).any(axis=0)]
        db = db[np.asarray([db['B(x)']!=elec[0], db['B(y)']!=elec[1]]).any(axis=0)]
        db = db[np.asarray([db['M(x)']!=elec[0], db['M(y)']!=elec[1]]).any(axis=0)]
        db = db[np.asarray([db['N(x)']!=elec[0], db['N(y)']!=elec[1]]).any(axis=0)]
    print('{} points removed (out of {})'.format(nbPoints-len(db), nbPoints))
    return db

def getLines(filename):
    headers = "Time	MeasID	DPID	Channel"
    endOfFile = "-----------------------"
    toFix = False
    file = open(filename,"r")
    for lineNb, line in enumerate(file):
        if headers in line:
            lineHeaders = lineNb
            if "		" in line:
                toFix = True
                print("The current file is broken. Trying to fix it")
        elif endOfFile in line:
            lineEnd = lineNb
    if toFix:
        fileInit = open(filename,"r").readlines()
        lineToEdit = fileInit[lineHeaders]
        while "		" in lineToEdit:
            lineToEdit = lineToEdit.replace("		", "	")
        fileInit[lineHeaders] = lineToEdit
        out = open(filename, "w")
        out.writelines(fileInit)
        out.close()

    return lineHeaders, lineEnd

def saveFiltered(filename, extension, db, headersLineInit, endOfFileInit):
    fileInit = open(filename,"r").readlines()
    filename = filename.replace(".txt","_"+extension+".txt")
    out = open(filename,"w", newline="")
    out.writelines(fileInit[:headersLineInit])
    out.writelines(db.to_csv(sep='\t',index=False))
    out.writelines(fileInit[endOfFileInit-1:])
    out.close()

if __name__ == "__main__":

    filenames = [f for f in listdir(directory) if isfile(join(directory, f))]

    for filename in filenames:
        lineHeaders, lineEnd = getLines(filename=join(directory,filename))
        DataBase = pd.read_csv(join(directory,filename), skiprows=lineHeaders, nrows=lineEnd-lineHeaders, sep='\t', index_col=False)
        DataBase = removeFaultyElectrode(DataBase, coordRemove)
        DataBase = DataBase[DataBase['Var(%)']<varMax]

        for i, yUsed in enumerate(yLines):
            DataBaseModif = extractLine(DataBase, y=yUsed)
            saveFiltered(join(directory,filename), nameLinesY[i], DataBaseModif, lineHeaders, lineEnd)

        for i, xUsed in enumerate(xLines):
            DataBaseModif = extractTransect(DataBase, x=xUsed)
            saveFiltered(join(directory,filename), nameLinesX[i], DataBaseModif, lineHeaders, lineEnd)