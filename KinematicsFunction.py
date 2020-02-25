import numpy as np  # Import NumPy
import math as m
import settings as settings

frequency = settings.frequency
omega = frequency*2*np.pi

def WingKinematics(offset, amplitude, time, depht):
    motion = offset + amplitude*m.sin(omega*time + depht)

    return motion