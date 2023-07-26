import numpy as np
from scipy.constants import c
path = r'C:\\Users\bengl\Documents\\PyHEADTAIL-playground\\IW2D-master\\examples\\outputs\\flat_chamber\\'
# Convert wake from IW2D to PyHT convention
A = np.loadtxt(path+r'ZxdipWLHC_1layersup_0layersdown4.00mm_test_TCDQ_Cu.dat', skiprows=1)
B = np.loadtxt(path+r'ZydipWLHC_1layersup_0layersdown4.00mm_test_TCDQ_Cu.dat', skiprows=1)
C = np.loadtxt(path+r'ZxquadWLHC_1layersup_0layersdown4.00mm_test_TCDQ_Cu.dat', skiprows=1)
D = np.loadtxt(path+r'ZyquadWLHC_1layersup_0layersdown4.00mm_test_TCDQ_Cu.dat', skiprows=1)
E = np.loadtxt(path+r'ZlongWLHC_1layersup_0layersdown4.00mm_test_TCDQ_Cu.dat', skiprows=1)

n_samples = A.shape[0]
W = np.zeros((n_samples, 6))

W[:,0] = A[:,0]
W[:,1] = A[:,1]
W[:,2] = B[:,1]
W[:,3] = C[:,1]
W[:,4] = D[:,1]
W[:,5] = E[:,1]

# PyHEADTAIL conversion
gamma = 32.9736
beta = np.sqrt( 1.-1./gamma**2)

# Time
W[:,0] /= (beta * c)
W[:,0] *= 1e9

# Transverse
W[:,1] /= 1e15
W[:,2] /= 1e15
W[:,3] /= 1e15
W[:,4] /= 1e15

# Longitudinal
W[:,5] /= 1e12

if np.isnan(W).any():
    print('Detected some NaNs, you should remove them')

# Save to new file
np.savetxt('MR_fullWake_resistiveWall_30GeV_PyHTConvention.dat', W)
