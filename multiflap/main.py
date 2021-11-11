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
from numpy.linalg import norm
from numpy import inf
from ms_package.lma_solver import Solver
from utilities.save_data import SaveData
from scipy import optimize
amplitude_hist = []
w_hist = []
def vertical_velocity(amplitude_shoulder):

# =============================================================================
# Generate the bird obj (setting and wing kinematics)
# =============================================================================

  #  bird_settings = SimulationsSettings(70.)

    mybird.shoulder.axis_z.amplitude = amplitude_shoulder

  #  bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
  #                     axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

  #  bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
  #                     axis_z=Joint(0.,0.,0.))

  #  mybird = BirdModel(shoulder=bird_shoulder,
  #                     elbow=bird_elbow,
  #                     wrist=bird_wrist,
  #                     settings=bird_settings)

  #  # Store the settings in a .txt file
    SaveData(sim_name).write_simulation_settings(mybird)


# =============================================================================
# Generate the multiple-shooting object
# =============================================================================

    x0 = initial_value[-1] # initial guess points

    ms_obj = MultipleShooting(x0, M = 2, period=0.25, t_steps=70, model = mybird)
# =============================================================================
# Solve the multiple-shooting eqs with LMA algorithm
# =============================================================================

    mysol = Solver(tolerance = 3e-5, ms_obj = ms_obj, max_iterations=25).lma()

    eigenvalues, eigenvectors = np.linalg.eig(mysol[4])
    sol_array = mysol[3].space
    sol_time = mysol[3].time
    u = sol_array[:,0]
    w = sol_array[:,1]
    theta = sol_array[:,3]
    initial_value.append(sol_array[0, :])
    initial_value.remove(initial_value[-2])
    # update initial condition
    u_fix = u*np.cos(theta) + w*np.sin(theta) # fixed frame velocity
    w_fix = u*np.sin(theta) - w*np.cos(theta)
    w_mean = np.mean(w_fix)

# =============================================================================
# Save the results in .npy format for the postprocessing
# =============================================================================
    ampli = amplitude_hist.append(amplitude_shoulder)
    amp_array = np.asarray(ampli)
    w_list = w_hist.append(w_mean)
    w_history = np.asarray(w_list)
    SaveData(sim_name).save_data('multipliers_LC', eigenvalues)
    SaveData(sim_name).save_data('solution_LC', sol_array)
    SaveData(sim_name).save_data('time_array_LC', sol_time)
    SaveData(sim_name).save_data('u_fix', u_fix)
    SaveData(sim_name).save_data('w_fix', w_fix)
    SaveData(sim_name).save_data('amp_shoulder', amp_array)
    SaveData(sim_name).save_data('w_history', w_history)
    print(w_mean)

    return w_mean

sweep_linspace = np.arange(22., 22.5, 0.5)
tail_linspace = np.arange(30., 30.5, 0.5)
initial_value = []
initial_value.append(np.array([18., 0.5, 0.1, 0.01]))
for i in range(len(sweep_linspace)):

    for j in range(len(tail_linspace)):
        sim_name = "t_"+str(tail_linspace[j])+"_s2_"+str(round(sweep_linspace[i], 2))
        SaveData(sim_name).make_folder()
        bird_settings = SimulationsSettings(tail_linspace[j])
        bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                                 axis_y=Joint(-np.deg2rad(sweep_linspace[i]),np.deg2rad(25),np.pi/2),
                                 axis_z=Joint(0,np.deg2rad(40),np.pi))


        bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                           axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

        bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                           axis_z=Joint(0.,0.,0.))

        mybird = BirdModel(shoulder=bird_shoulder,
                           elbow=bird_elbow,
                           wrist=bird_wrist,
                           settings=bird_settings)
        root = optimize.newton(vertical_velocity, np.deg2rad(34), tol = 5e-4)
        SaveData(sim_name).save_data('amplitude', root)
        #try:
        #    root = optimize.newton(vertical_velocity, np.deg2rad(34), tol = 5e-4)
        #    SaveData(sim_name).save_data('amplitude', root)
        #except:
        #    print("error")
            #SaveData(sim_name).remove_directory()
            #SaveData(sim_name).remove_contents()
            #SaveData(sim_name).remove_directory()
