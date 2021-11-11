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
from utilities.LimitCycleForces import RetrievingAero
from scipy import optimize
import scipy.interpolate as interp
count = 0


list_interpolation = np.load('./interpolation_values.npy')
def vertical_velocity(x):
    twist = x
    print("twist shoulder x: ", twist)
# =============================================================================
# Generate the bird obj (setting and wing kinematics)
# =============================================================================
    global count
    count = count + 1
    mybird.elbow.axis_x.offset = twist



    simulation_saving.write_simulation_settings(mybird)

# =============================================================================
# Generate the multiple-shooting object
# =============================================================================
    x0 = initial_value[-1]

    ms_obj = MultipleShooting(x0, M = 2, period=0.25, t_steps=90, model = mybird)
# =============================================================================
# Solve the multiple-shooting eqs with LMA algorithm
# =============================================================================

    mysol = Solver(tolerance = 1e-4, ms_obj = ms_obj, max_iterations=25).lma()

    eigenvalues, eigenvectors = np.linalg.eig(mysol[4])
    sol_array = mysol[3].space
    sol_time = mysol[3].time

    u = sol_array[:,0]
    w = sol_array[:,1]
    q = sol_array[:,2]
    theta = sol_array[:,3]
    initial_value.append(sol_array[0, :])
    #initial_value.remove(initial_value[-2])
    post_running = RetrievingAero(mybird, sol_time)


    [Fx,
    Fy,
    Fz,
    My,
    F_tail,
    M_wing,
    M_tail,
    M_drag,
    M_lift,
    P_tot,
    P_pro,
    aoa] =  post_running.postrun_aeroforces(sol_array)
    P_par = post_running.parasitic_power(sol_array)
    power_tot = np.mean(P_tot)
    # =============================================================================
    # Save the results in .npy format for the postprocessing
    # =============================================================================

    # update initial condition
    u_fix = u*np.cos(theta) + w*np.sin(theta) # fixed frame velocity
    w_fix = u*np.sin(theta) - w*np.cos(theta)
    u_mean = np.mean(u_fix)
    w_mean = np.mean(w_fix)

    simulation_saving.save_data('multipliers_LC', eigenvalues)
    simulation_saving.save_data('solution_LC', sol_array)
    simulation_saving.save_data('time_array_LC', sol_time)
    simulation_saving.save_data('u_fix', u_fix)
    simulation_saving.save_data('w_fix', w_fix)
    simulation_saving.save_data('lift', Fy)
    simulation_saving.save_data('wing_moment', M_wing)
    simulation_saving.save_data('tail_moment', M_tail)
    simulation_saving.save_data('angle_of_attack', aoa)
    simulation_saving.save_data('power', power_tot)
    print("vertical velocity: ", w_mean)
    print("power: ", power_tot)

    return w_mean


initial_value = []

tail_opening = 0.

sweep_offset = 10
bird_settings = SimulationsSettings(tail_opening)
simulation_saving = SaveData(project_name = 'test',
                             simulation_id = 'test_elbow')

simulation_saving.make_folder_sim()
simulation_saving.make_folder()
twist_x = 0.1
twist_elbow = 0.


x_int = list_interpolation[:, 0]
y_int = list_interpolation[:, 1]
u_int = list_interpolation[:, 2]
w_int = list_interpolation[:, 3]
q_int = list_interpolation[:, 4]
t_int = list_interpolation[:, 5]
o_int = list_interpolation[:, 6]

interpolator_u = interp.NearestNDInterpolator(list(zip(x_int, y_int)), u_int)
interpolator_w = interp.NearestNDInterpolator(list(zip(x_int, y_int)), w_int)
interpolator_q = interp.NearestNDInterpolator(list(zip(x_int, y_int)), q_int)
interpolator_t = interp.NearestNDInterpolator(list(zip(x_int, y_int)), t_int)
interpolator_o = interp.NearestNDInterpolator(list(zip(x_int, y_int)), o_int)

u_guess = interpolator_u(tail_opening, sweep_offset)
w_guess = interpolator_w(tail_opening, sweep_offset)
q_guess = interpolator_q(tail_opening, sweep_offset)
t_guess = interpolator_t(tail_opening, sweep_offset)
o_guess = interpolator_o(tail_opening, sweep_offset)


initial_value.append(np.array([u_guess, w_guess, q_guess, t_guess]))
bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                        axis_y=Joint(-np.deg2rad(sweep_offset),np.deg2rad(14),np.pi/2),
                         axis_z=Joint(0.2,np.deg2rad(37),np.pi))

bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                   axis_y=Joint(np.deg2rad(25),np.deg2rad(25),-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.deg2rad(30),np.deg2rad(30),np.pi/2),
                   axis_z=Joint(0.,0.,0.))

#bird_shoulder = Shoulder(axis_x=Joint(0.1, 0.014, -np.pi/2),
#                         axis_y=Joint(-np.deg2rad(24),np.deg2rad(20),np.pi/2),
#                         axis_z=Joint(0,np.deg2rad(40),np.pi))

#bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
#                   axis_y=Joint(np.deg2rad(30),np.deg2rad(30),-np.pi/2))
#
#bird_wrist = Wrist(axis_y=Joint(-np.deg2rad(30),np.deg2rad(30),np.pi/2),
#                   axis_z=Joint(0.,0.,0.))

mybird = BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist, settings=bird_settings, mass=2.2)

simulation_saving = SaveData(project_name = 'increase_mass', simulation_id = 'tail3_'+'{:.{prec}f}'.format(np.rad2deg(mybird.tail_opening), prec=3)+'_sweep'+'{:.{prec}f}'.format(np.rad2deg(mybird.shoulder.axis_y.offset), prec=3))
simulation_saving.make_folder_sim()
simulation_saving.make_folder()

root = optimize.newton(vertical_velocity, np.deg2rad(o_guess) ,tol = 5e-4)
