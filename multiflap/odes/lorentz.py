'''
file lorentz.py
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


class Lorentz:
    def __init__(self, sigma=None, r=None, b=None):

        self.sigma = sigma
        self.r = r
        self.b = b
        self.dimension = 3

    def dynamics(self, x0, t):

        """ODE system
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values
            t: time

        Outputs:
            x_dot: velocity vector
        """
        x, y, z = x0
        dxdt = self.sigma*(y-x)
        dydt = x*(self.r-z)-y
        dzdt = x*y - self.b*z
        vel_array = np.array([dxdt, dydt, dzdt], float)
        return vel_array

    def get_stability_matrix(self, x0, t):

        """
        Stability matrix of the ODE system

        Inputs:
            x0: initial condition
        Outputs:
            A: Stability matrix evaluated at x0. (dxd) dimension
            A[i, j] = dv[i]/dx[j]
        """
        x, y, z = x0
        A_matrix = np.array([[-self.sigma,  self.sigma, 0],
                            [(self.r-z), -1, -x],
                            [y, x, -self.b]], float)

        return A_matrix
