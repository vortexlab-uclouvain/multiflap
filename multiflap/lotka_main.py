'''
file rossler_main.py
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
from  odes.lotka import lotka
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting_period import MultipleShootingPeriod
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import SolverPeriod

x = [10., 10.]

time_array = np.linspace(0, 180, 90000)
mymodel = lotka(a=0.2, b=0.2, c=5.7)

ms_obj =  MultipleShootingPeriod(x, M=2, period_guess= 5., t_steps=50000, model=mymodel, option_jacobian='numerical')

mysol = SolverPeriod(ms_obj = ms_obj).lma()


eigenvalues, eigenvectors = np.linalg.eig(jac)



plt.legend()
plt.show()
