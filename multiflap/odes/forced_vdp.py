'''
file forced_vdp.py
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

"""
Example case adopted from:

    Practical Bifurcation and Stability Analysis, page 347
    Seydel R.
    Forced Van der Pol oscillator


    gamma = 0.5 ; lambda = 1.8
    p:q = 1:3 (from Seydel results). Meaning the period T=10.4719755

    This example does not use on purpose the stability matrix, in order to
    calculate the Jacoian numerically only.
"""

class ForcedVanDerPol:
    def __init__(self, lam=1.8 , gamma = 0.5, delta=4):
        self.lam = lam
        self.gamma = gamma
        self.delta = delta
        self.dimension=2
    def dynamics(self, x0, t):

        """ODE system
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values
            t: time

        Outputs:
            x_dot: velocity vector
        """
        y1, y2 = x0
        dy1_dt = y2
        dy2_dt = self.delta*(1 - y1**2)*y2 -y1 + self.gamma*(np.cos(self.lam*t))

        vel_array = np.array([dy1_dt, dy2_dt], float)
        return vel_array
