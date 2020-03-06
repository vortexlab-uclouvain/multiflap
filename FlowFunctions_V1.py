import numpy as np  # Import NumPy
import RungeKutta as rk
from FlappingForcesDev import FlappingForces
import settings as settings
import scipy.integrate as ode
import time
from math import sin, cos
dim = 4

f = settings.frequency

def birdEqn_py(t, ssp, **kinematics):
        """
        State space velocity function for the Equation of Motion of the longitudinal plane
    
        Inputs:
        ssp: State space vector
        ssp = (u, w, q, theta)
        t: Time
        
        Outputs:
        vel: Time derivative of ssp.
        """
        #Parameters:
        g = 9.81
        mass = 1.2
        #Read inputs:
        u, w, q, theta  = ssp  # Read state space points
        u = ssp[0]
        w = ssp[1]
        q = ssp[2]
        theta = ssp[3]
        # Longitudinal equation of motion:
        [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FlappingForces(t, u, w, q, theta, **kinematics)
        dudt = -q*w - g*sin(theta) - Fz/mass 
        dwdt = q*u + g*cos(theta) - Fy/mass - F_tail/mass
        dqdt =  My/0.1
        dthetadt = q                # Collect Equations of motion in a single NumPy array:
    
        vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector            
        return vel

def StabilityMatrix(t, ssp, **kinematics):
    """
    Stability matrix for the Rossler flow

    Inputs:
    ssp: State space vector. dx1 NumPy array: ssp = [x, y, z]
    Outputs:
    A: Stability matrix evaluated at ssp. dxd NumPy array
       A[i, j] = del Velocity[i] / del ssp[j]
    """
    g = 9.81
    m = 1.2
    perturbation = 1e-3
    u, w, q, theta = ssp
    [Fx, Fy, Fz, My, F_tail, _, _, _, _] = FlappingForces(t, u, w, q, theta, **kinematics)
    [Fxu, Fyu, Fzu, Myu, F_tailu, _, _, _, _] = FlappingForces(t, u + u*perturbation, w, q, theta, **kinematics)
    [Fxw, Fyw, Fzw, Myw, F_tailw, _, _, _, _] = FlappingForces(t, u ,w + w*perturbation, q, theta, **kinematics)
    [Fxq, Fyq, Fzq, Myq, F_tailq, _, _, _, _] = FlappingForces(t, u ,w , q + q*perturbation, theta, **kinematics)

    # Derivatives of Fz with respect to the state space variables
    dFzu_dU = (Fzu - Fz)/(u*perturbation)
    dFzw_dW = (Fzw - Fz)/(w*perturbation)
    dFzq_dq = (Fzq - Fz)/(q*perturbation)
    
    # Derivatives of Fy with respect to the state space variables
    dFyu_dU = (Fyu - Fy)/(u*perturbation)
    dFyw_dW = (Fyw - Fy)/(w*perturbation)
    dFyq_dq = (Fyq - Fy)/(q*perturbation)
    
    # Derivatives of F_tail with respect to the state space variables    
    dFytail_du = (F_tailu - F_tail)/(u*perturbation)
    dFytail_dw = (F_tailw - F_tail)/(w*perturbation)
    dFytail_dq = (F_tailq - F_tail)/(q*perturbation)
    
    # Derivatives of My with respect to the state space variables        
    dMy_du = (Myu - My)/(u*perturbation)
    dMy_dw = (Myw - My)/(w*perturbation)
    dMy_dq = (Myq - My)/(q*perturbation)

    A = np.array([[-dFzu_dU/m, -q - dFzw_dW/m, -w - dFzq_dq/m, -g*cos(theta)],
                  [q - dFyu_dU/m - dFytail_du/m, -dFyw_dW/m - dFytail_dw/m, u - dFyq_dq/m - dFytail_dq/m, -g*sin(theta)],
                  [dMy_du/0.1, dMy_dw/0.1, dMy_dq/0.1, 0],
                  [0, 0, 1, 0]], float) 
    return A

def JacobianVelocity(sspJacobian, t, **kinematics):
    """
    Velocity function for the Jacobian integration

    Inputs:
    sspJacobian: (d+d^2)x1 dimensional state space vector including both the
                 state space itself and the tangent space
    t: Time. Has no effect on the function, we have it as an input so that our
       ODE would be compatible for use with generic integrators from
       scipy.integrate

    Outputs:
    velJ = (d+d^2)x1 dimensional velocity vector
    """
    
    # Prendo il vettore sspJacobian (d + d^2) e la parte d la riempio, mentre 
    # la parte d^2 la trasformo nel Jacobiano (d x d)
    ssp = sspJacobian[0:dim]  # First three (or d in general) elements form the original state
                            # space vector
    J = sspJacobian[dim:].reshape((dim, dim))  # Last nine elements corresponds to
                                         # the elements of Jacobian.
    #We used numpy.reshape function to reshape d^2 dimensional vector which
    #hold the elements of Jacobian into a dxd matrix.
    #See
    #http://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
    #for the reference for numpy.reshape function

    velJ = np.zeros(np.size(sspJacobian))  # Initiate the velocity vector as a
                                           # vector of same size with
                                           # sspJacobian (d + d^2)
#    velJ[0:dim] = birdEqn_py(t, ssp)
    velJ[0:dim] = Velocity(ssp, t, **kinematics)

   
    #Last dxd elements of the velJ are determined by the action of
    #stability matrix on the current value of the Jacobian:

    velTangent = np.dot(StabilityMatrix(t, ssp, **kinematics), J)  # Velocity matrix (calculated in the ssp
                                                  # points)  for the tangent space

    velJ[dim:] = np.reshape(velTangent, dim**2)  # Another use of numpy.reshape, here
                                          # to convert from dxd to d^2

    return velJ


def Jacobian(ssp, t_initial, integration_time, **kinematics):
    """
    Jacobian function for the trajectory started on ssp, evolved for time t

    Inputs:
    ssp: Initial state space point. dx1 NumPy array: ssp = [x, y, z]
    t: Integration time
    Outputs:
    J: Jacobian of trajectory f^t(ssp). dxd NumPy array
    """
    #CONSTRUCT THIS FUNCTION
    #Hint: See the Jacobian calculation in CycleStability.py
    Jacobian0 = np.identity(dim) #Initial condition for Jacobian matrix
#    Jacobian0[0,1] = 0.1
    sspJacobian0 = np.zeros(dim + dim ** 2)  # Initiate
    sspJacobian0[0:dim] = ssp  # First 3 elemenets
    sspJacobian0[dim:] = np.reshape(Jacobian0, dim**2)  # Remaining 9 elements
    #print(sspJacobian0)
    t_Final = t_initial + integration_time  # Final time
    Nt = 15  # Number of time points to be used in the integration
    tArray = np.linspace(t_initial, t_Final, Nt)  # Time array for solution
    start_jac = time.time()
#    sspJacobianSolution = ode.solve_ivp(JacobianVelocity,[t_initial, t_Final], sspJacobian0, 'RK45')
    sspJacobianSolution = rk.RK4(JacobianVelocity, sspJacobian0, tArray, **kinematics)
    end_jac = time.time()
    print("Jacobian time ", (end_jac-start_jac))

#    sspJacobianSolution = sspJacobianSolution.y.T
    #Read the Jacobian for the periodic orbit:
    J = sspJacobianSolution[-1, dim:].reshape((dim, dim))
    return J

def Velocity(ssp, t, **kinematics):
        """
        State space velocity function for the Equation of Motion of the longitudinal plane

        Inputs:
        ssp: State space vector
        ssp = (u, w, q, theta)
        t: Time
        
        Outputs:
        vel: Time derivative of ssp.
        """
        #Parameters:
        g = 9.81
        mass = 1.2
        #Read inputs:
        u, w, q, theta  = ssp  # Read state space points
        u = ssp[0]
        w = ssp[1]
        q = ssp[2]
        theta = ssp[3]
        # Longitudinal equation of motion:
        [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FlappingForces(t, u, w, q, theta, **kinematics)
        dudt = -q*w - g*np.sin(theta) - Fz/mass 
        dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
        dqdt =  My/0.1
        dthetadt = q
        # Collect Equations of motion in a single NumPy array:

        vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector
        
        return vel


def Flow(ssp0, initial_time, deltat, time_steps, **kinematics):
    """
    Lagrangian description of the flow:
    This function integrates Rossler equation starting at ssp0 for deltat, and
    returns the final state space point.
    Inputs:
    ssp0: Initial state space point
    deltat: Integration time
    Outputs:
    sspdeltat: Final state space point
    """
    #Following numerical integration will return a 2 by 3(=d) solution array
    #where first row contains initial point ssp0, and the last row contains
    #final point
    tInitial =  initial_time            # Initial time
    tFinal = initial_time + deltat     # Final time
    Nt = time_steps                     # Number of time points to be used in the integration

    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array for solution

    sspSolution = rk.RK4(Velocity, ssp0, tArray, **kinematics) # RK 
#    sspSolution = ode.solve_ivp(birdEqn_py, [tInitial, tFinal], ssp0,'RK23', max_step = deltat/Nt)
#    sspSolution = (sspSolution.y).T
    sspdeltat = sspSolution[-1, :]  # Read the final point to sspdeltat
    return sspdeltat, sspSolution


def JacobianNumerical(ssp, initial_time, integration_time, **kinematics):
    
    """
    Finite difference evaluation of Jacobian
    Flow is perturbed in all directions of the phase space, and the generic
    component of the Jacobian is calculated by finite difference as follow:
        
    dF[i]/dx[j] = (F^t(x_perturbed)[i] - F^t(x)[i])/perturbation
    
    Jacobian = dF[i]/dx[j] (dim x dim) matrix
    
    epsilon = value of the perturbation   
    """
    time_steps = 50
    # -------------------------------------------------------------------------
    #  Initialization of the Jacobian Matrix
    # -------------------------------------------------------------------------
    
    Jacobian = np.zeros((dim,dim))
    
    # -------------------------------------------------------------------------
    # Set the numerical perturbation over the direction of the flow
    # -------------------------------------------------------------------------
    
    epsilon = 1e-3

    # ------------------------------------------------------------------------- 
    # Unperturbed initial condition
    # -------------------------------------------------------------------------
    
    ssp0 = ssp
    
    # -------------------------------------------------------------------------
    # Finite difference scheme for Jacobian evaluation
    # -------------------------------------------------------------------------
    
    print("... Running Jacobian Function")
    for j in range (dim):
        perturbation = np.zeros(dim)

        perturbation[j] = perturbation[j] + epsilon 
        ssp_pert = ssp + perturbation
        [vel, _] = Flow(ssp0, initial_time, integration_time, time_steps, **kinematics)
        [vel_pert, _] =  Flow(ssp_pert, initial_time, integration_time, time_steps, **kinematics)
        for i in range (dim):
            Jacobian[i,j] = (vel_pert[i] - vel[i])/perturbation[j]
        
        
    print("... Jacobian Calculated")
    return Jacobian


"""
The following block integrates numerically the equations of motion coupled with
the lifting line method. This is totally disconnected from the multi-shooting
code, it's just coupled numerical integration over a certain time.
"""    
if __name__ == "__main__":
    
    force_retrieving =  False
#    case_name = 'TestCase12b'

    if force_retrieving ==  True:    
        case_name = 'TestOffset'
        results_directory = '/Users/gducci/UCL/PROJECT/Simulations/LevelFlight/Storing_Results/Sweep_25/Opening_45'
        periodic_orbit_filename = results_directory+'/complete_solution.npy'
        periodic_orbit = np.load(periodic_orbit_filename)
        u_0 =  periodic_orbit[0,0][0]    # U-velocity initial condition
        w_0 =  periodic_orbit[0,0][1]         # W-velocity initial condition
        q_0 =  periodic_orbit[0,0][2]    # Q-velocity initial condition
        theta_0 = periodic_orbit[0,0][3]     # Theta-angle initial condition
        ssp0 = np.array([u_0+(0.1*u_0), w_0+(0.1*w_0), q_0, theta_0+(0.1*theta_0)], float) # Vector of initial conditions
        periodic_orbit = periodic_orbit.reshape(-1, periodic_orbit.shape[2])

    else:
        u_0 =  20.   # U-velocity initial condition
        w_0 = .0          # W-velocity initial condition
        q_0 = 0. # Q-velocity initial condition
        theta_0 = 0.  # Theta-angle initial condition
        ssp0 = np.array([u_0, w_0, q_0, theta_0], float) # Vector of initial conditions
        

    # -------------------------------------------------------------------------
    #  Definition of the time array
    # -------------------------------------------------------------------------
    period_number = 1       # Number of periods over which integrate the EoM
    tInitial = 0                      # Initial time
    tFinal = period_number*(1/f)      # Final time (This is one period)
    Nt = 30*period_number                   # Discretisation of time array
    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array
    
    sspSolution_V0 = rk.RK4(Velocity, ssp0,
                            tArray,
                            off_shoulder_y=-0.2, off_shoulder_x=0.5)
    
    import matplotlib.pyplot as plt
    plt.plot(tArray,  sspSolution_V0[:,1], '.', color = 'red', label = 'RK2')

###    plt.plot(tArray,  sspSolution_RK4[:,0], color = 'blue', label = 'RK4')
#    plt.plot(tArray,  sspSolution[:,1], color = 'green', label = 'RK4')
##
##
###    plt.plot(tArray,  USolution_RK3)
#    plt.show()
#    start_Jac = time.time()
#    Jacobian = JacobianNumerical(ssp0, tFinal)
#    end_Jac = time.time()
#    time_Jac = end_Jac - start_Jac
#    print(Jacobian)
#    print("RK2 = ", time_Jac)
#



