import numpy as np
from numpy import loadtxt
from scipy import interpolate
#
#global polar
#polar = np.load('polar.npy')
#
#def lift_coefficient(alpha):
#    if alpha < polar[-1,0]:
#        aoa = polar[:,0]
#        cl_polar = polar[:,1]
#        f = interpolate.interp1d(aoa[:], cl_polar[:], kind='linear')
#        cl = f(alpha)
#    else:
#        cl = polar[-1,1]
#    return cl