'''
file goodwin_main.py
@author Gianmarco Ducci, Marta del Olmo
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
from  odes.goodwin_model import GoodwinModel
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting_period import MultipleShootingPeriod
from matplotlib.ticker import FormatStrFormatter
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import SolverPeriod
import matplotlib as mpl

x = [0.1, 1, 1]

time_array = np.linspace(0, 250, 90000)
mymodel = GoodwinModel()

ms_obj =  MultipleShootingPeriod(x,
                                 M=40,
                                 period_guess= 23.,
                                 t_steps=200,
                                 model=mymodel,
                                 integrator='odeint', option_jacobian='numerical')

mysol = SolverPeriod(ms_obj = ms_obj).lma()

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)



sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

ssp_pert = odeint(GoodwinModel(perturbation=True).dynamics, sol_array[-1,:], time_array)
ssp_unpert = odeint(GoodwinModel(perturbation=False).dynamics, sol_array[-1,:], time_array)
pert = GoodwinModel().square(time_array, amp=2, t_0 = 150)

fig1 = plt.figure(figsize=(12,8))
ax1 = fig1.add_subplot(312)
ax1.plot(time_array, ssp_pert[:,0], label='x')
ax1.plot(time_array, ssp_pert[:,1], label='y')
ax1.plot(time_array, ssp_pert[:,2], label='z')
ax1.legend()
ax1.set_xlabel('time (h)'); ax1.set_ylabel('Goodwin variables')
ax2 = fig1.add_subplot(311)
ax2.plot(time_array, pert)
ax2.set_xlabel('time (h)'); ax2.set_ylabel('Gaussian perturbation')
ax3 = fig1.add_subplot(313)
ax3.plot(time_array, ssp_pert[:,0] - ssp_unpert[:,0], label='x_pert - x_unpert')
ax3.plot(time_array, ssp_pert[:,1] - ssp_unpert[:,1], label='y_pert - y_unpert')
ax3.plot(time_array, ssp_pert[:,2] - ssp_unpert[:,2], label='z_pert - z_unpert')
ax3.legend()
ax3.set_xlabel('time (h)'); ax3.set_ylabel('pert-unpert')
plt.show()
