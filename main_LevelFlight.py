import numpy as np
import FlowFunctions_V1 as func
import MultiShooting_newscheme as multi2
from RungeKutta import RK4
import settings as settings
from LimitCycleForces import ForceRetrieving
from scipy import optimize
import os
from scipy import interpolate

# =============================================================================
# Step 0: Create the Main folder where all results will be stored in
# =============================================================================

path_directory = "/Users/gducci/UCL/PROJECT/Simulations/ResultsTail/LevelSimulations/"          # Type here the path where you wish to create your folder

simulation_name = input("Type here the name of the folder: ") # Input from terminal

simulation_directory = path_directory+simulation_name
try:
    os.mkdir(simulation_directory) # Create target Directory
except FileExistsError:
    print("... Directory already exists, please type another name")
tolerance_MultiShooting = float(input("Type here the tolerance value: "))     # Input name from terminal
tolerance_RootFinding = 5e-4 # Tolerance value for finding the mean vertical velocity
initial_guessed_value=[]
initial_guessed_value.append(np.array([16.0, 0.5, 0.1, 0.01]))

def vertical_velocity(amplitude_shoulder, sweep_amplitude, sweep_offset):
    print("Starting Newton scheme with amplitude shoulder z: ", np.rad2deg(amplitude_shoulder))
#    global COUNT
#    COUNT += 1

    dim = func.dim
    period = 1/(settings.frequency)
    M = multi2.M
    states_stack = np.zeros([M,dim])
    tau = (period)/(M-1)
    op_tail = 0.    # deg
    # Initial points distribution
    
    states_stack[0,0:] = initial_guessed_value[-1]
    
    for i in range (1,M):
        [states_stack[i,0:], _] = func.Flow(states_stack[i-1,:], i*tau, tau, 50, 
                                            amp_shoulder_y=sweep_amplitude,
                                            amp_shoulder_z=amplitude_shoulder,
                                            off_shoulder_y=sweep_offset,
                                            tail_opening=np.deg2rad(op_tail))
    
    # Call the multiple shooting routine
    
    xfin, ptlist, error, complete_solution, Jacobian_semigroup = multi2.MultiShootingScheme(states_stack, 
                                                                                            tau, 
                                                                                            100, 
                                                                                            tolerance_MultiShooting, 
                                                                                            simulation_directory,
                                                                                            amp_shoulder_y=sweep_amplitude,
                                                                                            amp_shoulder_z=amplitude_shoulder,
                                                                                            off_shoulder_y=sweep_offset,
                                                                                            tail_opening=np.deg2rad(op_tail))
    
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
    print("Mean Vertical Velocity = ", mean_vertical_velocity)
    
    opening_folder=simulation_directory+"/SweepAmplitude_"+str(int(round(abs(np.rad2deg(sweep_amplitude)))))+"/SweepOff_Neg"+str(int(round(np.rad2deg(sweep_offset))))
    try:  
        os.mkdir(opening_folder)
    except OSError:  
        print ("Creation of the directory %s failed" % opening_folder)
    else:  
        print ("Successfully created the directory %s " % opening_folder)

    LevelFlight = []
    LevelFlight.append(['Amplitude shoulder z', amplitude_shoulder])
    LevelFlight.append(['Amplitude sweep angle', sweep_amplitude])
    LevelFlight.append(['Offset sweep angle', sweep_offset])
    LevelFlight.append(['Mean Forward Flight velocity', np.mean(u_hor)])
    ptlist=np.asarray(ptlist)
# =============================================================================
#         Saving results
# =============================================================================
    np.save(opening_folder+'/LevelFlight_Result', LevelFlight)    
    np.save(opening_folder+'/Error_History', np.asanyarray(error))    
    np.save(opening_folder+'/outfile_ptlist', ptlist)
    np.save(opening_folder+'/complete_solution', complete_solution)
    np.save(opening_folder+'/final_points', complete_solution)
    np.save(opening_folder+'/initial_points', states_stack)
    eigenValues_SG, eigenVectors_SG = np.linalg.eig(Jacobian_semigroup)
    np.save(opening_folder+'/outfile_JacobianEigenvalues_SG', eigenValues_SG)
    np.save(opening_folder+'/outfile_JacobianEigenvector_SG', eigenVectors_SG)
    eigenValues, eigenVectors = np.linalg.eig(Jacobian_semigroup)
    np.save(opening_folder+'/outfile_JacobianEigenvalues', eigenValues)
    np.save(opening_folder+'/outfile_JacobianEigenvector', eigenVectors)
    print("...Retrieving Aerodynamic forces and moments")
    [Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, 
     Moment_tail, Moment_drag, Moment_lift] = ForceRetrieving(complete_solution, 
                                                             amp_shoulder_y=sweep_amplitude,
                                                             amp_shoulder_z=amplitude_shoulder,
                                                             off_shoulder_y=sweep_offset,
                                                             tail_opening=np.deg2rad(op_tail))
    
    np.save(opening_folder+'/Lift_coupled_v2', Fy)
    np.save(opening_folder+'/Drag_coupled_v2', Fz)
    np.save(opening_folder+'/Force_tail', F_tail)
    np.save(opening_folder+'/Moment_total', Moment_total)
    np.save(opening_folder+'/Moment_wing', Moment_wing)
    np.save(opening_folder+'/Moment_lift', Moment_lift)
    np.save(opening_folder+'/Moment_drag', Moment_drag)
    np.save(opening_folder+'/Moment_tail', Moment_tail)
    np.save(opening_folder+'/Amplitude_shoulder_z', np.rad2deg(amplitude_shoulder))
    print('Saving results completed')
# =============================================================================
#       Writing DataFile.txt
# =============================================================================
    simulation_info = open(opening_folder+'/DataFile.txt', 'w+')
    simulation_info.write('Period: '+str(period)+'\n')
    simulation_info.write('FRAMES COORDINATES (with respect to the aerodynamic reference system)'+'\n')
    simulation_info.write('Wing frame position: '+str(settings.wingframe_position)+'\n')
    simulation_info.write('Tail frame position: '+str(settings.wingframe_position_tail)+'\n')
    simulation_info.write('Tail opening identically zero (no tail config) '+'\n')
    simulation_info.write('... Multiple shooting settings ...'+'\n')
    simulation_info.write('Number of points in shooting: '+str(M)+'\n')
    simulation_info.write('Number of Iterations: '+str(np.size(error))+'\n')
    simulation_info.write('Error: '+str((error[-1]))+'\n')
    simulation_info.write('Floquet Multipliers: '+str(eigenValues)+'\n')
    simulation_info.write('|$Lambda_1$|: '+str(np.linalg.norm(eigenValues[0]))+'\n')
    simulation_info.write('|$Lambda_2$|: '+str(np.linalg.norm(eigenValues[1]))+'\n')
    simulation_info.write('|$Lambda_3$|: '+str(np.linalg.norm(eigenValues[2]))+'\n')
    simulation_info.write('|$Lambda_4$|: '+str(np.linalg.norm(eigenValues[3]))+'\n')
    simulation_info.write('Average vertical velocity: '+str(mean_vertical_velocity)+'\n')
    simulation_info.write('Average forward flight velocity: '+str(np.mean(u_hor))+'\n')
    simulation_info.write('Shoulder amp z: '+str(np.rad2deg(amplitude_shoulder))+'\n')

    simulation_info.close()
    print('DataFile correctly written')
    return mean_vertical_velocity

# =============================================================================
# Interpolator to get close to the amplitude
# =============================================================================
    

#def amplitude_interpolator(sweep_angle, opening_angle):
#    points = np.load('/Users/gducci/UCL/PROJECT/AmplitudeInterpolation/points.npy')
#    amplitudes = np.load('/Users/gducci/UCL/PROJECT/AmplitudeInterpolation/amplitudes.npy')
#    f = interpolate.interp2d(points[:,0], points[:,1], amplitudes, kind='linear')
#    amplitude = f(sweep_angle, opening_angle)
#    return amplitude
#
#test=amplitude_interpolator(30,40.5)
# =============================================================================
# Define the range of tail opening over which iterate
# =============================================================================
sweep_offset_min = 19
sweep_offset_max = 19
points = int((sweep_offset_max-sweep_offset_min)//1)
sweep_offset_range = np.linspace(sweep_offset_min, sweep_offset_max, (points+1), endpoint=True)

sweep_min = 20
sweep_max = 20
points = int((sweep_max-sweep_min)//5)
sweep_amplitude_range = np.linspace(sweep_min, sweep_max, points+1, endpoint=True)
points = int((50-sweep_min)//5)

amplitude_range = np.linspace(15, 50, points+1, endpoint=True)


for i in range(len(sweep_amplitude_range)):
    sweep_folder=simulation_directory+"/SweepAmplitude_"+str(int(sweep_amplitude_range[i]))
    try:  
        os.mkdir(sweep_folder)
    except OSError:  
        print ("Creation of the directory %s failed" % sweep_folder)
    else:  
        print ("Successfully created the directory %s " % sweep_folder)
    COUNT =0
    for j in range(len(sweep_offset_range)):       
#        amplitude_shoulder = float(amplitude_interpolator(sweep_range[i],tail_range[j]))
        root = optimize.newton(vertical_velocity, np.deg2rad(42), 
                               args=(np.deg2rad(sweep_amplitude_range[i]),np.deg2rad(-sweep_offset_range[j]),),
                               tol=tolerance_RootFinding)
#w=[]
#for i in range(len(amplitude_range)):
#    w_mean = vertical_velocity(np.deg2rad(amplitude_range[i]), sweep_min, opening_min)
#    w.append(w_mean)
#np.save(simulation_name+'/W_MEAN', w_mean)
#np.save(simulation_name+'/amplitudes', amplitude_range)

print("Simulation done")