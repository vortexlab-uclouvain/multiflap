import numpy as np
import FlowFunctions_V1 as func
import MultiShooting_newscheme as multi2
from RungeKutta import RK4, RK2
from settings import SimulationsSettings
from LimitCycleForces import ForceRetrieving
import os
import DataFile_writing as data_writing
import SaveFile as saving
from bird import BirdLiftingline, Joint
# -------------------------------------------------------------------------
# Creation of the Folder where simulation data will be stored
# It is simply asked to type from terminal the name the folder you wish to create
# -------------------------------------------------------------------------

path_directory = "../Simulations/setting_class/"          # Type here the path where you wish to create your folder

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
settings = SimulationsSettings()

bc = BirdLiftingline()

bc.shoulder.axis_x = Joint(0.2, 0.014, -np.pi/2)
bc.shoulder.axis_y  = Joint(-np.deg2rad(19), np.deg2rad(20), np.pi/2)
bc.shoulder.axis_z  = Joint(0., np.deg2rad(42), np.pi)

bc.elbow.axis_x = Joint(0., np.pi/6, -np.pi/2)
bc.elbow.axis_y = Joint(np.pi/6, np.pi/6, -np.pi/2)

bc.wrist.axis_y = Joint(-np.pi/6, np.pi/6, np.pi/2)
bc.wrist.axis_z = Joint(0., 0., 0.)

dim = func.dim
period = 1/(settings.frequency)
M = multi2.M
states_stack = np.zeros([M,dim])
tau = (period)/(M-1)

# Type here the first guessed point, that will be used to calculate the other (M-1) points

states_stack[0,0:] =[16., 0.5, 0.1, 0.01 ]

# Automatic routine to extract the remaining M-1 points for the flow. 
# Note this is not always the best way to guess points
tail_op = np.deg2rad(0)
settings.amplitude_shoulder_y = np.deg2rad(20)
for i in range (1,M):
    [states_stack[i,0:], _] = func.Flow(states_stack[i-1,:], (i-1)*tau, tau, 10)

guessed_points = np.copy(states_stack)


data_writing.DataFileParameters(simulation_directory)

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
                                                                                        simulation_directory)
results_directory = simulation_directory+"/Results_2"
os.mkdir(results_directory) # Create target Directory

ptlist=np.asarray(ptlist)


tArray = np.linspace(0, period, len(complete_solution))

"""
Printing floquet multipliers
"""

eigenValues_SG, eigenVectors_SG = np.linalg.eig(Jacobian_semigroup)
data_writing.DataFileResults(simulation_directory, error, eigenValues_SG)

[_,trajectory_recomputed] = func.Flow(xfin[0], 0, 0.25, multi2.time_steps*(multi2.M - 1))

print('Calculating Jacobian Numerical')
Jacobian_numerical = func.JacobianNumerical(xfin[0, :], 0, period)

eigenValues, eigenVectors = np.linalg.eig(Jacobian_numerical)
print("...Retrieving Aerodynamic forces and moments")
[Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift] = ForceRetrieving(complete_solution)

saving.SaveData(results_directory,ptlist, complete_solution,
             states_stack,
             error,
             eigenValues_SG,
             eigenVectors_SG,
             eigenValues,
             eigenVectors,
             trajectory_recomputed,
             Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift)

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
M
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
