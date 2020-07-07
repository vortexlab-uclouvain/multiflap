import numpy as np
import FlowFunctions_V1 as func
import MultiShooting_newscheme as multi2
from RungeKutta import RK4
import settings as settings
from LimitCycleForces import ForceRetrieving
from scipy import optimize
import os
from scipy import interpolate
import DataFile_writing as data_writing
import SaveFile as saving

# ==============================================================
# Step 0: Create the Main folder where all results will be stored in
# ===============================================================

path_directory = "/Users/gducci/UCL/PROJECT/Simulations/TailMap"


simulation_directory = path_directory

tolerance_MultiShooting = float(input("Type here the tolerance value: "))
tolerance_RootFinding = 5e-4  # Tolerance value mean vertical velocity
initial_guessed_value = []
initial_guessed_value.append(np.array([18.2, -1.9, -0.1, -0.115]))


def vertical_velocity(amplitude_shoulder, sweep_amplitude, sweep_offset):
    print("Amplitude shoulder z: ", np.rad2deg(amplitude_shoulder))
#    global COUNT
#    COUNT += 1

    dim = func.dim
    period = 1/(settings.frequency)
    M = multi2.M
    states_stack = np.zeros([M, dim])
    tau = (period)/(M-1)
    op_tail = 60.    # deg
    # Initial points distribution
    opening_folder = simulation_directory+"/TailOpening_"+str(round(op_tail))
    try:
        os.mkdir(opening_folder)
    except OSError:
        print("Creation of the directory %s failed" % opening_folder)
    else:
        print("Successfully created the directory %s " % opening_folder)
    settings.amplitude_shoulder_z = amplitude_shoulder
    settings.amplitude_shoulder_y = sweep_amplitude
    settings.offset_shoulder_y = sweep_offset
    settings.tail_opening = np.deg2rad(op_tail)
    states_stack[0, 0:] = initial_guessed_value[-1]

    for i in range(1, M):
        [states_stack[i, 0:], _] = func.Flow(states_stack[i-1, :], (i-1)*tau, tau, 50)

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
    print("Mean Vertical Velocity = ", mean_vertical_velocity)

    storing_results = opening_folder+"/SweepOffset_"+str(round(np.rad2deg(sweep_offset)))
    try:
        os.mkdir(storing_results)
    except OSError:
        print("Creation of the directory %s failed" % storing_results)
    else:
        print("Successfully created the directory %s " % storing_results)

    LevelFlight = []
    LevelFlight.append(['Amplitude shoulder z', amplitude_shoulder])
    LevelFlight.append(['Amplitude sweep angle', sweep_amplitude])
    LevelFlight.append(['Offset sweep angle', sweep_offset])
    LevelFlight.append(['Mean Forward Flight velocity', np.mean(u_hor)])
    ptlist = np.asarray(ptlist)
# ================================================================
# Saving results
# ================================================================
    eigenValues_SG, eigenVectors_SG = np.linalg.eig(Jacobian_semigroup)
    eigenValues, eigenVectors = np.linalg.eig(Jacobian_semigroup)
    print("...Retrieving Aerodynamic forces and moments")
    [Fx, Fy, Fz, Moment_total, F_tail, Moment_wing,
     Moment_tail, Moment_drag, Moment_lift] = ForceRetrieving(complete_solution)

    saving.SaveData(storing_results, ptlist, complete_solution,
                     states_stack,
                     error,
                     eigenValues_SG,
                     eigenVectors_SG,
                     eigenValues,
                     eigenVectors,
                     Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift)

    print('Saving results completed')
# ================================================================
#       Writing DataFile.txt
# ================================================================
    data_writing.DataFileParameters(storing_results)
    data_writing.DataFileResults_Level(storing_results, error, eigenValues, mean_vertical_velocity, u_hor, amplitude_shoulder)

    print('DataFile correctly written')
    return mean_vertical_velocity

# ================================================================
# Interpolator to get close to the amplitude
# ================================================================


#def amplitude_interpolator(sweep_angle, opening_angle):
#    points = np.load('/Users/gducci/UCL/PROJECT/AmplitudeInterpolation/points.npy')
#    amplitudes = np.load('/Users/gducci/UCL/PROJECT/AmplitudeInterpolation/amplitudes.npy')
#    f = interpolate.interp2d(points[:,0], points[:,1], amplitudes, kind='linear')
#    amplitude = f(sweep_angle, opening_angle)
#    return amplitude
#
#test=amplitude_interpolator(30,40.5)
# ================================================================
# Define the range of tail opening over which iterate
# ================================================================
sweep_offset_min = 22
sweep_offset_max = 26
points = int((sweep_offset_max-sweep_offset_min)//1)
sweep_offset_range = np.linspace(sweep_offset_min, sweep_offset_max, (points+1), endpoint=True)

sweep_min = 20
sweep_max = 20
points = int((sweep_max-sweep_min)//5)
sweep_amplitude_range = np.linspace(sweep_min, sweep_max, points+1, endpoint=True)
points = int((50-sweep_min)//5)

amplitude_range = np.linspace(15, 50, points+1, endpoint=True)


for i in range(len(sweep_amplitude_range)):
    COUNT =0
    for j in range(len(sweep_offset_range)):
#        amplitude_shoulder = float(amplitude_interpolator(sweep_range[i],tail_range[j]))
        root = optimize.newton(vertical_velocity, np.deg2rad(42),
                               args=(np.deg2rad(sweep_amplitude_range[i]),np.deg2rad(-sweep_offset_range[j]),),
                               tol=tolerance_RootFinding)

print("Simulation done")
