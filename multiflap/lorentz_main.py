'''
file lorentz_main.py
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
from  odes.lorentz import MyModel
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting_period import MultipleShooting
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import Solver
from mpl_toolkits.mplot3d import Axes3D

x = [10., 10., 3.6]

mymodel = MyModel(a=10, b=28, c=8/3)

ms_obj =  MultipleShooting(x, M=2, period_guess= 10., t_steps=50000, model=mymodel)

odes = mymodel.dynamics

t_array = np.linspace(0., 50, 50000)
sol = rk4(odes, x, t_array)

mysol = Solver(ms_obj = ms_obj).lma(5., '/Users/gducci/UCL/PROJECT/Simulations/class_test')

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)

sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

sol_integration = sol.x

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
ax.plot(sol_integration[:, 0],
        sol_integration[:, 1],
        sol_integration[:, 2], alpha = 0.3, color='red')

ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'b')
plt.show()
