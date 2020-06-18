import numpy as np
import FlowFunctions_V1 as func
import MultiShooting_newscheme as multi2
import settings as settings
from LimitCycleForces import ForceRetrieving
from scipy import optimize
import os
from scipy import interpolate
import DataFile_writing as data_writing
import SaveFile as saving

# =============================================================================
# Step 0: Create the Main folder where all results will be stored in
# =============================================================================

path_directory = "/Users/gducci/UCL/PROJECT/Simulations/ResultsTail/LevelSimulations/"

simulation_name = input("Type here the name of the folder: ")  # Input terminal

simulation_directory = path_directory+simulation_name
try:
    os.mkdir(simulation_directory)  # Create target Directory
except FileExistsError:
    print("... Directory already exists, please type another name")
tolerance_MultiShooting = float(input("Type here the tolerance value: "))
tolerance_RootFinding = 5e-4  # Tolerance value mean vertical velocity
initial_guessed_value = []
initial_guessed_value.append(np.array([18.2, -1.9, -0.1, -0.115]))
COUNT = 0


def vertical_velocity(amplitude_shoulder, sweep_offset):
    print("Starting Newton scheme with amplitude shoulder z: ", np.rad2deg(amplitude_shoulder))
    global COUNT

    COUNT += 1

    dim = func.dim
    period = 1/(settings.frequency)
    M = multi2.M
    states_stack = np.zeros([M, dim])
    tau = (period)/(M-1)
    op_tail = 20.    # deg
    # Initial points distribution
    settings.amplitude_shoulder_z = amplitude_shoulder
    settings.offset_shoulder_y = sweep_offset
    settings.tail_opening = np.deg2rad(op_tail)
    states_stack[0, 0:] = initial_guessed_value[-1]

    for i in range(1, M):
        [states_stack[i, 0:], _] = func.Flow(states_stack[i-1, :],
                                             (i-1)*tau, tau, 50)

    # Call the multiple shooting routine

    xfin, ptlist, error, complete_solution, Jacobian_semigroup = multi2.MultiShootingScheme(states_stack, 
                                                                                            tau, 
                                                                                            100, 
                                                                                            tolerance_MultiShooting, 
                                                                                            simulation_directory)

    periodic_orbit = complete_solution.reshape(-1, complete_solution.shape[2])
    initial_guessed_value.append(periodic_orbit[0])
    initial_guessed_value.remove(initial_guessed_value[-2])
    velocity_U = periodic_orbit[:, 0]
    velocity_W = periodic_orbit[:, 1]
    Q = periodic_orbit[:, 2]
    Theta = periodic_orbit[:, 3]
    u_hor = velocity_U*np.cos(Theta) + velocity_W*np.sin(Theta)
    v_vert = velocity_U*np.sin(Theta) - velocity_W*np.cos(Theta)
    mean_vertical_velocity = np.mean(v_vert)
    mean_horizontal_velocity = np.mean(u_hor)
    print("Mean Vertical Velocity = ", mean_vertical_velocity)

    opening_folder = simulation_directory+"/Results_"+str(COUNT)
    try:
        os.mkdir(opening_folder)
    except OSError:
        print("Creation of the directory %s failed" % opening_folder)
    else:
        print("Successfully created the directory %s " % opening_folder)

    LevelFlight = []
    LevelFlight.append(['Amplitude shoulder z', amplitude_shoulder])
    LevelFlight.append(['Offset sweep angle', sweep_offset])
    LevelFlight.append(['Mean Forward Flight velocity', np.mean(u_hor)])
    ptlist = np.asarray(ptlist)
# ===========================================================
# Saving results
# ===========================================================
    eigenValues_SG, eigenVectors_SG = np.linalg.eig(Jacobian_semigroup)
    eigenValues, eigenVectors = np.linalg.eig(Jacobian_semigroup)
    print("...Retrieving Aerodynamic forces and moments")
    [Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail,
     Moment_drag, Moment_lift] = ForceRetrieving(complete_solution)

    saving.SaveData(opening_folder, ptlist, complete_solution,
                    states_stack,
                    error,
                    eigenValues_SG,
                    eigenVectors_SG,
                    eigenValues,
                    eigenVectors,
                    Fx, Fy, Fz, Moment_total, F_tail,
                    Moment_wing, Moment_tail, Moment_drag, Moment_lift)

    print('Saving results completed')
# ===========================================================
# Writing Datafile.txt
# ===========================================================
    data_writing.DataFileParameters(opening_folder)
    data_writing.DataFileResults_Level(opening_folder, error, eigenValues,
                                       mean_vertical_velocity, u_hor,
                                       amplitude_shoulder)

    print('DataFile correctly written')
    return mean_horizontal_velocity, mean_vertical_velocity


amplitude = np.linspace(36.8, 37.2, 5)

offset_y = np.linspace(-21.1, -20.8, 4)

fixed_frame_velocities = []

kinematics = []

for i in range(len(amplitude)):

    for j in range(len(offset_y)):

        u_ff, w_ff = vertical_velocity(np.deg2rad(amplitude[i]), np.deg2rad(offset_y[j]))

        sol = [u_ff, w_ff]

        fixed_frame_velocities.append(sol)

        kinematics = [amplitude[i], offset_y[j]]

        np.save(simulation_directory+"/Results_"+str(COUNT)+"/sol", sol)

        np.save(simulation_directory+"/Results_"+str(COUNT)+"/kinematics", kinematics)
