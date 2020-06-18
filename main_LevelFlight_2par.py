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

# =============================================================================
# Step 0: Create the Main folder where all results will be stored in
# =============================================================================

path_directory = "/Users/gducci/UCL/PROJECT/Simulations/2ParametersOptimization/"          # Type here the path where you wish to create your folder

simulation_name = input("Type here the name of the folder: ") # Input from terminal

simulation_directory = path_directory+simulation_name
try:
    os.mkdir(simulation_directory) # Create target Directory
except FileExistsError:
    print("... Directory already exists, please type another name")
tolerance_MultiShooting = float(input("Type here the tolerance value: "))     # Input name from terminal
tolerance_RootFinding = 5e-4 # Tolerance value for finding the mean vertical velocity
initial_guessed_value=[]
initial_guessed_value.append(np.array([18.2, -1.9, -0.1, -0.115]))

U_ForFlight = 16
def my_fun(params):
#    print("Starting Newton scheme with amplitude shoulder z: ", np.rad2deg(amplitude_shoulder))
#    global COUNT
#    COUNT += 1
    print(np.tan(params))
    params_unpacked = np.tan(params)
    dim = func.dim
    period = 1/(settings.frequency)
    M = multi2.M
    states_stack = np.zeros([M,dim])
    tau = (period)/(M-1)
    op_tail = 20.    # deg
    # Initial points distribution
    settings.amplitude_shoulder_z = params_unpacked[0]
    settings.amplitude_shoulder_y = np.deg2rad(20)
    settings.offset_shoulder_y = params_unpacked[1]
    settings.tail_opening = np.deg2rad(op_tail)
    states_stack[0,0:] = initial_guessed_value[-1]
    
    for i in range (1,M):
        [states_stack[i,0:], _] = func.Flow(states_stack[i-1,:], (i-1)*tau, tau, 50)
    
    # Call the multiple shooting routine
    
    xfin, ptlist, error, complete_solution, Jacobian_semigroup = multi2.MultiShootingScheme(states_stack, 
                                                                                            tau, 
                                                                                            100, 
                                                                                            tolerance_MultiShooting, 
                                                                                            simulation_directory)
    
    periodic_orbit = complete_solution.reshape(-1, complete_solution.shape[2])
    initial_guessed_value.append(periodic_orbit[0])
    initial_guessed_value.remove(initial_guessed_value[-2])
    velocity_U = periodic_orbit[:,0]
    velocity_W = periodic_orbit[:,1]
    Q = periodic_orbit[:,2]
    Theta = periodic_orbit[:,3]
    u_hor = velocity_U*np.cos(Theta) + velocity_W*np.sin(Theta)
    v_vert = velocity_U*np.sin(Theta) - velocity_W*np.cos(Theta)
    mean_vertical_velocity = np.mean(v_vert)
    mean_horizontal_velocity = np.mean(u_hor)
    print("Mean Vertical Velocity = ", mean_vertical_velocity)
    print("Mean Horizontal Velocity = ", mean_horizontal_velocity)
    
    opening_folder=simulation_directory+"/Results"
    try:  
        os.mkdir(opening_folder)
    except OSError:  
        print ("Creation of the directory %s failed" % opening_folder)
    else:  
        print ("Successfully created the directory %s " % opening_folder)

    LevelFlight = []
    LevelFlight.append(['Amplitude shoulder z', params[0]])
#    LevelFlight.append(['Amplitude sweep angle', sweep_amplitude])
    LevelFlight.append(['Offset sweep angle', params[1]])
    LevelFlight.append(['Mean Forward Flight velocity', np.mean(u_hor)])
    ptlist=np.asarray(ptlist)
# =============================================================================
#         Saving results
# =============================================================================
    eigenValues_SG, eigenVectors_SG = np.linalg.eig(Jacobian_semigroup)
    eigenValues, eigenVectors = np.linalg.eig(Jacobian_semigroup)
    print("...Retrieving Aerodynamic forces and moments")
    [Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, 
     Moment_tail, Moment_drag, Moment_lift] = ForceRetrieving(complete_solution)
    
    saving.SaveData(opening_folder,ptlist, complete_solution,
             states_stack,
             error,
             eigenValues_SG,
             eigenVectors_SG,
             eigenValues,
             eigenVectors,
             Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift)

    print('Saving results completed')
# =============================================================================
#       Writing DataFile.txt
# =============================================================================
    data_writing.DataFileParameters(opening_folder)
    data_writing.DataFileResults_Level(opening_folder, error, eigenValues, mean_vertical_velocity, u_hor, params[0])
    
    print('DataFile correctly written')
    velocity_norm = np.sqrt(((mean_horizontal_velocity/U_ForFlight - 1)**2) + ((mean_vertical_velocity/(1.2*9.81*0.25))**2))
    print("norm = ", velocity_norm)
    return velocity_norm

amp_shoulder_z = np.arctan(np.deg2rad(37))

sweep_offset_min = np.arctan(-np.deg2rad(20.9))

xyf0 = np.array([amp_shoulder_z, sweep_offset_min])
xyf_bounds=((np.arctan(np.deg2rad(36)), np.arctan(np.deg2rad(38))),( np.arctan(-np.deg2rad(22)), np.arctan(-np.deg2rad(18))))

root = optimize.basinhopping(my_fun, xyf0, niter=100, T=1.0, stepsize=0.5)

print("Simulation done")