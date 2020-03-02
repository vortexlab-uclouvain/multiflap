import numpy as np
import FlowFunctions_V1 as func
import MultiShooting_newscheme as multi2
from RungeKutta import RK4
import settings as settings
from LimitCycleForces import ForceRetrieving
import os

# -------------------------------------------------------------------------
# Creation of the Folder where simulation data will be stored
# It is simply asked to type from terminal the name the folder you wish to create
# -------------------------------------------------------------------------

path_directory = "../Simulations/NoIteration/"          # Type here the path where you wish to create your folder

simulation_name = input("Type here the name of the folder: ")      # Input name from terminal
simulation_directory = path_directory+simulation_name

try:
    os.mkdir(simulation_directory) # Create target Directory
    print("... Directory Created ") 
except FileExistsError:
    print("... Directory already exists, please type another name")

tolerance = float(input("Type here the tolerance value: "))     # Input name from terminal

"""
Parameters of the flow:
    dim = dimension of the ODE system, this is read from the FlowFunction file    
    period = Period of the orbit. This is in the scheme used, known a priori    
    M = number of the guessed points for the multiple-shooting scheme
    states_stack = initialization of the matrix that stores the coordinates of the guessed points (M x dim)
    tau = time interval between two consecutives points
"""

dim = func.dim
period = 1/(settings.frequency)
M = multi2.M
states_stack = np.zeros([M,dim])
tau = (period)/(M-1)

# Type here the first guessed point, that will be used to calculate the other (M-1) points

states_stack[0,0:] =[20.9, -2.8, -0.09, -0.17]

# Automatic routine to extract the remaining M-1 points for the flow. 
# Note this is not always the best way to guess points
amplitude_shoulder = np.deg2rad(50)
sweep = np.deg2rad(25)
offset_shoulder_y = -np.deg2rad(29)
tail_op = np.deg2rad(0)

for i in range (1,M):
    [states_stack[i,0:], _] = func.Flow(states_stack[i-1,:], i*tau, tau, 50, amp_shoulder_y=sweep,
                                        amp_shoulder_z=amplitude_shoulder,
                                        off_shoulder_y=offset_shoulder_y,
                                        tail_opening=tail_op)
 
# Keep the guessed points in memory

guessed_points = np.copy(states_stack)

# =============================================================================
# Writing txt info simulation
# =============================================================================
simulation_info = open(simulation_directory+'/DataFile.txt', 'w+')
simulation_info.write('Period: '+str(period)+'\n')
simulation_info.write('FRAMES COORDINATES (with respect to the aerodynamic reference system)'+'\n')
simulation_info.write('Wing frame position: '+str(settings.wingframe_position)+'\n')
simulation_info.write('Tail frame position: '+str(settings.wingframe_position_tail)+'\n')
simulation_info.write('Tail opening: '+str(np.rad2deg(tail_op))+'\n')
simulation_info.write('Tail chord: '+str(settings.tail_length)+'\n')
simulation_info.write('... Multiple shooting settings ...'+'\n')
simulation_info.write('Number of points in shooting: '+str(M)+'\n')
simulation_info.write('Initial Guess: '+str(states_stack[0,0:])+'\n')

"""
MultiShootingScheme

---------------------------------------------------------------------------
---------------------------------------------------------------------------    

    Args: 
           states_stack: Coordinates of the points at each time step     
           tau:          Time interval between two consecutives points   
           maxIter:      Maximum number of iteration, defined in the main
           tol:          Tolerance before exiting the loop              
          
---------------------------------------------------------------------------
---------------------------------------------------------------------------    
          
    Return:
           x_fin:        Coordinates of the points that belong to the Limit cycle. 
           ptlist:       List of the point coordinates for post-processing
           error:        List of the error for every step of the iteration
          
---------------------------------------------------------------------------
---------------------------------------------------------------------------       
"""
xfin, ptlist, error, complete_solution, Jacobian_semigroup = multi2.MultiShootingScheme(states_stack, 
                                                                                        tau, 
                                                                                        100, 
                                                                                        tolerance, 
                                                                                        simulation_directory,
                                                                                        amp_shoulder_y=sweep,
                                                                                        amp_shoulder_z=amplitude_shoulder,
                                                                                        off_shoulder_y=offset_shoulder_y,
                                                                                        tail_opening=tail_op)


results_directory = simulation_directory+"/Results"
os.mkdir(results_directory) # Create target Directory

ptlist=np.asarray(ptlist)
np.save(results_directory+'/outfile_ptlist', ptlist)
np.save(results_directory+'/complete_solution', complete_solution)
np.save(results_directory+'/final_points', complete_solution)
np.save(results_directory+'/initial_points', states_stack)


tArray = np.linspace(0, period, len(complete_solution))

"""
Printing floquet multipliers
"""

eigenValues_SG, eigenVectors_SG = np.linalg.eig(Jacobian_semigroup)
np.save(results_directory+'/outfile_JacobianEigenvalues_SG', eigenValues_SG)
np.save(results_directory+'/outfile_JacobianEigenvector_SG', eigenVectors_SG)

#Jacobian = func.JacobianNumerical(xfin[0, :], 0, period)

eigenValues, eigenVectors = np.linalg.eig(Jacobian_semigroup)
np.save(results_directory+'/outfile_JacobianEigenvalues', eigenValues)
np.save(results_directory+'/outfile_JacobianEigenvector', eigenVectors)
print("...Retrieving Aerodynamic forces and moments")
[Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift] = ForceRetrieving(complete_solution,
                                                                                        amp_shoulder_y=sweep,
                                                                                        amp_shoulder_z=amplitude_shoulder,
                                                                                        off_shoulder_y=offset_shoulder_y,
                                                                                        tail_opening=tail_op)
np.save(results_directory+'/Lift_coupled_v2', Fy)
np.save(results_directory+'/Drag_coupled_v2', Fz)
np.save(results_directory+'/Force_tail', F_tail)
np.save(results_directory+'/Moment_total', Moment_total)
np.save(results_directory+'/Moment_wing', Moment_wing)
np.save(results_directory+'/Moment_lift', Moment_lift)
np.save(results_directory+'/Moment_drag', Moment_drag)
np.save(results_directory+'/Moment_tail', Moment_tail)

simulation_info.write('Number of Iterations: '+str(np.size(error))+'\n')
simulation_info.write('Error: '+str((error[-1]))+'\n')
simulation_info.write('Floquet Multipliers: '+str(eigenValues)+'\n')
simulation_info.write('|$Lambda_1$|: '+str(np.linalg.norm(eigenValues[0]))+'\n')
simulation_info.write('|$Lambda_2$|: '+str(np.linalg.norm(eigenValues[1]))+'\n')
simulation_info.write('|$Lambda_3$|: '+str(np.linalg.norm(eigenValues[2]))+'\n')
simulation_info.write('|$Lambda_4$|: '+str(np.linalg.norm(eigenValues[3]))+'\n')
simulation_info.write('Shoulder z amplitude: '+str(np.rad2deg(amplitude_shoulder))+'\n')
simulation_info.close() 

print("Eigenvalues are, ", eigenValues)
"""
Plotting the final solution
"""
trajectory = complete_solution.reshape(-1, complete_solution.shape[2])
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker

# Create a figure instance
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

# Set axis label
ax.set_xlabel('$x$')  
ax.set_ylabel('$y$')  
ax.set_zlabel('$z$') 


# Plot the periodic orbit:
ax.plot(trajectory[:, 0],
        trajectory[:, 1],
        trajectory[:, 2],color = 'b')

# Plot the final M-points
ax.scatter(xfin[:,0],
           xfin[:,1],
           xfin[:,2], color='red', label='Final points')

# Plot the initial guessed M-points
ax.scatter(states_stack[:,0], 
           states_stack[:,1], 
           states_stack[:,2], color = 'green', label='Guessed points')
plt.legend(loc='upper left', numpoints = 1 )

fig2 = plt.figure(2)
ax2 = fig2.gca()
ax2.set_aspect('equal')
ax2.set_xlabel('$Re$', fontsize=18)  
ax2.set_ylabel('$Im$', fontsize=18)  
ax2.xaxis.set_tick_params(labelsize=10)

circle = np.linspace(0,2*np.pi,101)
ax2.plot(np.cos(circle),np.sin(circle))
plt.grid()
ax2.plot(eigenValues.real, eigenValues.imag,'ro' )
ax2.plot(eigenValues_SG.real, eigenValues_SG.imag,'x' )
plt.show()
