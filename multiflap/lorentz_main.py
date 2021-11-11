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
from  odes.lorentz import Lorentz
from ms_package.multiple_shooting_period import MultipleShootingPeriod
from ms_package.lma_solver_period import SolverPeriod
from ms_package.lyapunov_exponents import LyapunovExponents
import matplotlib.pyplot as plt
x = [-5, -0., 20]

# Creating the object containing the eqs
mymodel = Lorentz(a=10, b=28, c=8/3)

#ms_obj =  MultipleShootingPeriod(x, M=10, period_guess= 8., t_steps=50000, model=mymodel, option_jacobian='analytical', integrator='odeint')
#
#mysol = SolverPeriod(ms_obj = ms_obj, tolerance=1e-4).lma()
#jac = mysol[4]
#
#sol_array = mysol[3].space
#sol_time = mysol[3].time
#period = sol_time[-1]
#
#x_new = [sol_array[0,0], sol_array[0,1], sol_array[0,2]]
ms_obj =  MultipleShootingPeriod(x, M=25, period_guess= 1.4, t_steps=200000, model=mymodel, option_jacobian='analytical', integrator='odeint')

mysol = SolverPeriod(ms_obj = ms_obj, tolerance=1e-4).lma()

sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]


t_f = 100
n = 10000

#time_array = np.linspace(20, 100, 2)

#hist = []
#for i in range (len(time_array)):
        ## calling the LyapunovExponents class
        #lambda_t = LyapunovExponents(x, n, time_array[i], mymodel).get_lyapunov_exponent()
        #hist.append(lambda_t)
        #print(lambda_t)

fig = plt.figure(figsize=(12, 5))
fig.suptitle("Period detected:  {}".format(period))
ax = fig.add_subplot(1, 2, 1, projection= '3d')
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'k')

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot( sol_time, sol_array[:,0], label = "x")
ax2.plot( sol_time, sol_array[:,1], label = "y")
ax2.plot( sol_time, sol_array[:,2], label = "z")

ax2.set_xlabel('$t$')
ax2.set_ylabel('$f(t)$')

plt.legend()
plt.savefig('/Users/gducci/Desktop/fig2.pdf')
plt.show()
print(period)
