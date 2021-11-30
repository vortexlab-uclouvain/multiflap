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
from ms_package.lyapunov_exponents import LyapunovExponents
from ms_package.multiple_shooting_period import MultipleShootingPeriod
from ms_package.lma_solver_period import SolverPeriod
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

#x = [15., 6.78, 3.6]
#x = [-5., -6.78, 3.6]

r_values = [24.73, 23, 22, 21]

fig1 = plt.figure(1)

for i in range(len(r_values)):
    # Creating the object containing the eqs
    mymodel = Lorentz(sigma=10, r=r_values[i], b=8/3)

    x = np.sqrt(mymodel.b*(mymodel.r-1))
    y = x
    z = mymodel.r-1

    fixed_points = np.array([x, y, z])

    initial_value = fixed_points*1.7
    t_f = 100
    n = 10000

    time_array = np.linspace(0, 2.5, 2000)

    #x_ms = [-7., -7., 3.6]
    integrator_sol = odeint(mymodel.dynamics, initial_value, time_array)
    r_hopf = mymodel.sigma*(mymodel.sigma + mymodel.b +3)/(mymodel.sigma - mymodel.b- 1)
    # %%
    ms_obj = MultipleShootingPeriod(initial_value, M=10, period_guess=1., t_steps=500, model=mymodel)
    mysol = SolverPeriod(ms_obj = ms_obj, max_iterations=200, tolerance=1e-8).lma()

    fixed_points = fsolve(mymodel.dynamics, [10, 5, 10], 0)
    sol_array = mysol[3].space
    sol_time = mysol[3].time

    print(r_hopf)
    x_p = sol_array[0, :] - sol_array[0,:]*0.1

    x_inner = sol_array[0, :] + sol_array[0,:]*0.5

    time_inner = np.linspace(0, 20, 2000)
    pert_sol = odeint(mymodel.dynamics, x_p, time_array)
    inner_sol = odeint(mymodel.dynamics, x_inner, time_inner)
    jac = mysol[4]

    eigenvalues, eigenvectors = np.linalg.eig(jac)
    # x_2 = fixed_points + 1
    #stable_int = odeint(mymodel.dynamics, x_2, time_array)

    # %%
    ax = fig1.gca(projection='3d')
    ax.plot(sol_array[:,0], sol_array[:,1], sol_array[:,2], label="r= "+str(mymodel.r))
    #ax.plot_trisurf(pert_sol[:,0], pert_sol[:,1], pert_sol[:,2], antialiased=True, alpha=0.3)
    #ax.plot(integrator_sol[:, 0], integrator_sol[:,1], integrator_sol[:,2])
    #ax.plot(pert_sol[:,0], pert_sol[:,1], pert_sol[:,2])
    #ax.plot(inner_sol[:,0], inner_sol[:,1], inner_sol[:,2])
    ax.scatter(-x, -y, z)
    #ax.scatter(x_2[0], x_2[1], x_2[2])
    ax.legend()
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
plt.savefig('/Users/gducci/Desktop/saddle_cycles.pdf', format = 'pdf')
plt.show()
#hist = []
#for i in range (len(time_array)):
#        # calling the LyapunovExponents class
#        lambda_t = LyapunovExponents(x, n, time_array[i], mymodel).get_lyapunov_exponent()
#        hist.append(lambda_t)
#        print(lambda_t)
#
#
#print(lambda_t)
