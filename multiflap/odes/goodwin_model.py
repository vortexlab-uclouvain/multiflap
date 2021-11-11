'''
file redox_oscillation.py
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

    A Robust Model for Circadian Redox Oscillations,
    del Olmo, M.; Kramer, A.; Herzel, H., Int. J. Mol. Sci., 2019, 20, 2368
    https://www.mdpi.com/1422-0067/20/9/2368

"""
class GoodwinModel:
    def __init__(self, k1=0.7, k2=0.45, k3=0.7, k4=0.35, k5=0.7, k6=0.35, 
                 K1=1, K2=1, K4=1, K6=1, n=7):

        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6 = k6
        self.K1 = K1
        self.K2 = K2
        self.K4 = K4
        self.K6 = K6
        self.n = n
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
        dx_dt = self.k1*(self.K1**self.n/(self.K1**self.n + z**self.n)) - \
                self.k2*x/(self.K2 + x)
        dy_dt = self.k3*x - self.k4*y/(self.K4 + y)
        dz_dt = self.k5*y - self.k6*z/(self.K6 + z)

        vel_array = np.array([dx_dt, dy_dt, dz_dt], float)
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
        A_matrix = np.array([[(-self.k2*(self.K2 + x) + self.k2*x)/(self.K2 + x)**2, 0., (-self.k1*(self.K1**self.n)*self.n*z**(self.n-1))/(self.K1**self.n + z**self.n)**2],
                      [self.k3,  (-self.k4*(self.K4 + y) + self.k4*y)/(self.K4 + y)**2, 0],
                      [0, self.k5, (-self.k6*(self.K6 + z) + self.k6*z)/(self.K6 + z)**2]], float)

        return A_matrix
