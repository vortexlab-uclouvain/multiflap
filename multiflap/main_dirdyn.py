'''
file main.py
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
from  ms_package.multiple_shooting import MultipleShooting
from aero_package.bird_model import BirdModel
from aero_package.settings import SimulationsSettings
from aero_package.bird import Shoulder, Elbow, Wrist, Joint
import time
from utilities.save_data import SaveData
from numpy.linalg import norm
from numpy import inf
from ms_package.lma_solver import Solver
import matplotlib.pyplot as plt

bird_settings = SimulationsSettings(0.)
sim_name = "paper_longrun"
SaveData(sim_name).make_folder()
# generate bird kinematics by calling "bird" module
bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                         axis_y=Joint(-np.deg2rad(19),np.deg2rad(20),np.pi/2),
                         axis_z=Joint(0,np.deg2rad(42),np.pi))
bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                   axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                   axis_z=Joint(0.,0.,0.))

mybird = BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist, settings=bird_settings)

sol = np.load("../results/paper_sim/paper_run_1/solution_LC.npy")
eigenvectors = np.load("../results/paper_sim/paper_run_1/eigenvectors.npy")
eigenvalues = np.load("../results/paper_sim/paper_run_1/floquet_multipliers.npy")
x = sol[0,:]
x_p = x - eigenvectors[:,-1]*0.1
# set initial guess for multiple-shooting scheme
periods = 5
time_steps = 60*periods
ms_obj = MultipleShooting(x, M = 2, period=0.25, t_steps=time_steps, model = mybird)
ms_obj_p = MultipleShooting(x_p, M = 2, period=0.25, t_steps=time_steps, model = mybird)

[_, initial_traj] = ms_obj.get_mappedpoint(x, 0, 0.25*periods)
[_, pert_traj] = ms_obj.get_mappedpoint(x_p, 0, 0.25*periods)

time_array_LC = initial_traj.t
solution_LC = initial_traj.x
time_array_pert = pert_traj.t
solution_pert = pert_traj.x

SaveData(sim_name).save_data('solution_unperturbed_s.npy', solution_LC)
SaveData(sim_name).save_data('time_unperturbed_s.npy', time_array_LC)
SaveData(sim_name).save_data('time_perturbed_s.npy', time_array_pert)
SaveData(sim_name).save_data('solution_perturbed_s.npy', solution_pert)
#fig1 = plt.figure(figsize=(12,5))
#ax1 = fig1.add_subplot(221)
#ax1.plot(time_array_LC, solution_LC[:,0], c='k')
#ax1.plot(time_array_pert, solution_pert[:,0], c='red')
#ax1.set_ylabel('$U$')
#ax1 = fig1.add_subplot(222)
#ax1.plot(time_array_LC, solution_LC[:,1], c='k')
#ax1.plot(time_array_pert, solution_pert[:,1], c='red')
#ax1.set_ylabel('$w$')
#ax1 = fig1.add_subplot(223)
#ax1.plot(time_array_LC, solution_LC[:,2], c='k')
#ax1.plot(time_array_pert, solution_pert[:,2], c='red')
#ax1.set_ylabel('$q$')
#ax1 = fig1.add_subplot(224)
#ax1.plot(time_array_LC, solution_LC[:,3], c='k')
#ax1.plot(time_array_pert, solution_pert[:,3], c='red')
#ax1.set_ylabel('$\Theta$')
##ax1.axhline(w_fix_mean, color='red', ls='--')
#plt.show()
#
