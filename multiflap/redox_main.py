'''
file redox_main.py
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
from  odes.redox_oscillation import RedoxModel
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting_period import MultipleShootingPeriod
from matplotlib.ticker import FormatStrFormatter
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import SolverPeriod
import matplotlib as mpl

x = [0.5, 0.5, 0.6, 0.2]

time_array = np.linspace(0, 180, 90000)
mymodel = RedoxModel()

ms_obj =  MultipleShootingPeriod(x, M=2, period_guess= 23., t_steps=50000, model=mymodel, integrator='odeint')

mysol = SolverPeriod(ms_obj = ms_obj).lma()

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)


sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

fig1 = plt.figure()
ax1 = fig1.gca()
ax1.plot( sol_time, sol_array[:,0], label = "D1")
ax1.plot( sol_time, sol_array[:,1], label = "D2")
ax1.plot( sol_time, sol_array[:,2], label = "R")
ax1.plot( sol_time, sol_array[:,3], label = "A")
plt.legend()
plt.savefig('../img/' + 'redox_oscillation.png', format = 'png' )
plt.show()

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.set_aspect('equal')
ax2.set_xlabel('$Re$', fontsize=20)  
ax2.set_ylabel('$Im$', fontsize=20)  
ax2.xaxis.set_tick_params(labelsize=10)
ax2.set_xlim(-1.1, 2.)
ax2.set_ylim(-1.2, 1.2)
ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1.))
ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1.))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.tick_params(labelsize=12)
circle = np.linspace(0,2*np.pi,101)
ax2.plot(np.cos(circle),np.sin(circle), linewidth=3.)
plt.grid(False)
ax2.scatter(eigenvalues.real, eigenvalues.imag , s=70,marker='.',color='red',facecolors='red', edgecolors='red', linewidth=1.8)
plt.gcf().subplots_adjust(left=0.16)
plt.savefig('../img/' + 'redox_multipliers.png', format = 'png' )
plt.show()
