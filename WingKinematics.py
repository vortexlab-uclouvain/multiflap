import numpy as np  # Import NumPy
import math as m

class WingKinematics:

    frequency = 4
    omega = frequency*2*np.pi
    def __init__(self, joint, axis, offset, amplitude, time, depht):
        self.joint = joint
        self.axis = axis
        self.offset = offset
        self.amplitude = amplitude
        self.time = time
        self.depht = depht
        self.motion = offset + amplitude*m.sin(WingKinematics.omega*time + depht)