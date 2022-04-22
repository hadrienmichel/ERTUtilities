import numpy as np
from scipy.spatial.transform import Rotation as rot
from matplotlib import pyplot as plt

def unitVector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angleBetween(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unitVector(v1)
    v2_u = unitVector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

MNT = np.loadtxt('AgatheData\PP9_ChessionAmont_MNT.xyz',delimiter=' ')
Electrodes = np.loadtxt('AgatheData\PP9_Electrodes.csv',skiprows=1,delimiter=',', usecols=(1,2,3))

# Reverse the electrodes (optional):
Electrodes = np.flipud(Electrodes)

MNT = MNT[MNT[:,-1] > 0, :] # Remove all unused points in MNT (marker = -9999)

# Define the transform:
# 1) Translation
Electrodes[:,:-1] -= Electrodes[0,:-1]
MNT[:,:-1] -= Electrodes[0,:-1]
# 2) Rotation
# Compute distance between the first and the last point (in XY plane):
dist = np.linalg.norm(Electrodes[0,:-1]-Electrodes[-1,:-1])
angle = angleBetween(Electrodes[-1,:-1], [dist, 0])
rotation = rot.from_rotvec([0, 0, -angle])
Electrodes = rotation.apply(Electrodes)
MNT = rotation.apply(MNT)
# transformMatrix = np.linalg.lstsq(Electrodes[(0,-1),:-1], np.asarray([[0, 0],[dist, 0]]))
# Electrodes[:,:-1] = np.matmul(Electrodes[:,:-1], transformMatrix[0])
# MNT[:,:-1] = np.matmul(MNT[:,:-1], transformMatrix[0])
plt.plot(Electrodes[:,0], Electrodes[:,1])
plt.show()

## Saving the electrodes positions and the MNT:
np.savetxt('AgatheData\PP9_ChessionAmont_MNT_Rotated.xyz', MNT, delimiter='\t')
np.savetxt('AgatheData\PP9_Electrodes_Rotated.xyz', Electrodes, delimiter='\t')
