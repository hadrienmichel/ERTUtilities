import numpy as np

## Files directory:
directory = './VSL_Data/Save_Line1/'
fileName = 'invResults.vtk'
fileNameChanges = 'invResultsChanged.vtk'

## Changes in positions:
# Points coordinates in the inverted profile
pointsInit = [[0,   0], # X Y position of the origin (0, 0)
              [126, 0]] # X Y position of the end of the profile
pointsInit = np.asarray(pointsInit)

# Points coordinates in real world
pointsFinal = [[204515.856,	93219.782],
               [204504.901,	93094.259]]
pointsFinal = np.array(pointsFinal)

zShift = 0 # If the topography is relative, you can add to get to the actual altitude

diffInit = pointsFinal[0,:]

changeYZ = True # Change the y to the z axis (for 2D profiles)

transformMatrix = np.linalg.lstsq(pointsInit, np.subtract(pointsFinal, diffInit))

print(np.matmul(pointsInit[0,:],transformMatrix[0]))
## Read and write *.vtk file:
fileInit = open(directory + fileName, 'r')
fileFinal = open(directory + fileNameChanges, 'w')

positions = False
for line in fileInit:
    if (positions) and (line.startswith('CELLS')):
        positions = False
    # Change the positions in the file:
    if positions:
        posCurr = np.asarray([float(s) for s in line.split('\t')])
        if changeYZ:
            posCurr[1], posCurr[2] = posCurr[2], posCurr[1]
        posFinal = posCurr
        posFinal[:2] = diffInit + np.matmul(posCurr[:2], transformMatrix[0])
        posFinal[2] += zShift
        fileFinal.writelines(f'{posFinal[0]}\t{posFinal[1]}\t{posFinal[2]}\n')
    else:
        fileFinal.writelines(line)
    if not(positions) and (line.startswith('POINTS')):
        positions = True

fileFinal.close()
fileInit.close()