import numpy as np
import math as m

frequency  = 4
omega = frequency*2*np.pi

class Joint:
    omega = frequency*2*np.pi

    def __init__(self, offset=None, amplitude=None, phase=None):
        self.offset = offset
        self.amplitude = amplitude
        self.phase = phase

    def motion_joint(self, t):

        motion = self.offset + self.amplitude*np.sin(omega*t + self.phase)
        return motion

class Shoulder():
    def __init__(self, axis_x=None, axis_y=None, axis_z=None):
        self.axis_x  = Joint()
        self.axis_y  = Joint()
        self.axis_z  = Joint()
        if isinstance(axis_x, Joint):
            self.axis_x = axis_x
        if isinstance(axis_y, Joint):
            self.axis_y = axis_y
        if isinstance(axis_z, Joint):
            self.axis_z = axis_z

class Elbow():
    def __init__(self, axis_x=None, axis_y=None):
        self.axis_x  = Joint()
        self.axis_y  = Joint()
        if isinstance(axis_x, Joint):
            self.axis_x = axis_x
        if isinstance(axis_y, Joint):
            self.axis_y = axis_y

class Wrist():
    def __init__(self, axis_y=None, axis_z=None):
        self.axis_y  = Joint()
        self.axis_z  = Joint()
        if isinstance(axis_y, Joint):
            self.axis_y = axis_y
        if isinstance(axis_z, Joint):
            self.axis_z = axis_z
