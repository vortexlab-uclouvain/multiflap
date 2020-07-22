'''
file bird.py
@author Gianmarco Ducci
@copyright Copyright © UCLouvain 2020

multiflap is a Python tool for finding periodic orbits and assess their stability via the Floquet multipliers.

Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>

List of the contributors to the development of multiflap, Description and complete License: see LICENSE and NOTICE files.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
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
