import numpy as np
from numpy import loadtxt
from scipy import interpolate

polar = np.load('polar.npy')

def lift_coefficient(alpha):
    if alpha < polar[-1,1]:
        aoa = polar[:,0]
        cl_polar = polar[:,1]
        f = interpolate.interp1d(aoa[:], cl_polar[:], kind='cubic')
        cl = f(alpha)
    else:
        cl = polar[-1,1]
    return cl


print(lift_coefficient(16.))