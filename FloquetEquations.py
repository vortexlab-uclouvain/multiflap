import numpy as np  # Import numpy
from FlappingForcesDev import FlappingForces 
import RungeKutta as rk
from FlowFunctions_V1 import JacobianNumerical
import scipy.integrate as ode
import time
a = None
c = None
dim = 5

def Velocity(ssp, t):
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
        u, w, q, theta, plus_one = ssp  # Read state space points
        u = ssp[0]
        w = ssp[1]
        q = ssp[2]
        theta = ssp[3]
        plus_one = ssp[4]
        # Longitudinal equation of motion:
        [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FlappingForces(t, u, w, q, theta)
        dudt = -q*w - g*np.sin(theta) - Fz/mass
        dwdt = q*u + g*np.cos(theta) - Fy/mass  - F_tail/mass
        dqdt =  My/0.1  
        dthetadt = q
        plus_onedt = 1
        # Collect Equations of motion in a single NumPy array:

        vel = np.array([dudt, dwdt, dqdt, dthetadt, plus_onedt], float)  # Velocity vector
        
        return vel

def StabilityMatrix(ssp, t):
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
    u, w, q, theta ,plus_one = ssp
    
    [Fx, Fy, Fz, My, F_tail, _, _, _, _] = FlappingForces(t, u, w, q, theta)
    [Fxu, Fyu, Fzu, Myu, F_tailu, _, _, _, _] = FlappingForces(t, u + u*perturbation, w, q, theta)
    [Fxw, Fyw, Fzw, Myw, F_tailw, _, _, _, _] = FlappingForces(t, u ,w + w*perturbation, q, theta)

    dFyu_dU = (Fyu - Fy)/(u*perturbation)
    dFzu_dU = (Fzu - Fz)/(u*perturbation)
    dFyw_dW = (Fyw - Fy)/(w*perturbation)
    dFzw_dW = (Fzw - Fz)/(w*perturbation)
    dFytail_du = (F_tailu - F_tail)/(u*perturbation)
    dFytail_dw = (F_tailw - F_tail)/(w*perturbation)
    dMy_du = (Myu - My)/(u*perturbation)
    dMy_dw = (Myw - My)/(w*perturbation)

    A = np.array([[-dFzu_dU/m, -q - dFzw_dW/m, -w, -g*np.cos(theta), 1],
                  [q - dFyu_dU/m - dFytail_du/m, -dFyw_dW/m - dFytail_dw/m, u, -g*np.sin(theta), 1],
                  [dMy_du/0.1, dMy_dw/0.1, 0, 0, 1],
                  [0, 0, 1, 0, 1],
                  [0 ,0,0,0,0]], float) 
    return A

def JacobianVelocity(t, sspJacobian):
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
    velJ[0:dim] = Velocity(ssp, t)
    
   
    #Last dxd elements of the velJ are determined by the action of
    #stability matrix on the current value of the Jacobian:

    velTangent = np.dot(StabilityMatrix(ssp, t), J)  # Velocity matrix (calculated in the ssp
                                                  # points)  for the tangent space

    velJ[dim:] = np.reshape(velTangent, dim**2)  # Another use of numpy.reshape, here
                                          # to convert from dxd to d^2

    return velJ


def Jacobian(ssp, t):
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
    tInitial = 0.  # Initial time
    tFinal = tInitial+t  # Final time
    Nt = 1000  # Number of time points to be used in the integration
    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array for solution

    sspJacobianSolution = ode.solve_ivp(JacobianVelocity,[tInitial, tFinal], sspJacobian0, 'RK23', max_step=tFinal/50)
    sspJacobianSolution = sspJacobianSolution.y.T
#    sspJacobianSolution = rk.RK2(JacobianVelocity, sspJacobian0, tArray)
    #Read the Jacobian for the periodic orbit:
    J = sspJacobianSolution[-1, dim:].reshape((dim, dim))
    return J


if __name__ == "__main__":
    #This block will be evaluated if this script is called as the main routine
    #and will be ignored if this file is imported from another script.
    #
    #This is a handy structure in Python which lets us test the functions
    #in a package

    #In order to test our integration routine, we are going to define Harmonic
    #Oscillator equations in a 2D state space:
    # -------------------------------------------------------------------------
    # ● Definition of the time array
    # -------------------------------------------------------------------------
    f = 4
    period_number = 1            # Number of periods over which integrate the EoM
    tInitial = 0.            # Initial time
    tFinal = tInitial+period_number*(1/f)      # Final time (This is one period)
    Nt = 50                      # Discretisation of time array
    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array

    # -------------------------------------------------------------------------
    # ● Initial conditions for numerical integration
    # -------------------------------------------------------------------------
			

    case_name = 'JacobianAnalytical'
    results_directory = '/Users/gducci/UCL/PROJECT/Simulations/JacobianTest/'+ case_name+'/Results'
    periodic_orbit_filename = results_directory+'/complete_solution.npy'
    periodic_orbit = np.load(periodic_orbit_filename)
    u_0 =  periodic_orbit[0,0][0]    # U-velocity initial condition
    w_0 =  periodic_orbit[0,0][1]         # W-velocity initial condition
    q_0 =  periodic_orbit[0,0][2]    # Q-velocity initial condition
    theta_0 = periodic_orbit[0,0][3]     # Theta-angle initial condition
    k = 1# Theta-angle initial condition
    ssp0 = np.array([u_0 , w_0, q_0, theta_0, k], float) # Vector of initial conditions
    # -------------------------------------------------------------------------
    # ● Call the numerical integration function
    # -------------------------------------------------------------------------
    start = time.time()

#    start_RK2 = time.time()
#    sspSolution_RK3 = ode.solve_ivp(Velocity,[tInitial, tFinal], ssp0, 'RK23', max_step=tFinal/50)
#    solution = sspSolution_RK3.y.T
    solution = rk.RK4(Velocity, ssp0, tArray)
    print(solution[0,0] - solution[-1,0])
#    end_RK2 = time.time()
#    time_RK2 = end_RK2 - start_RK2
#    print("RK2 = ", time_RK2)
#    jacobian = Jacobian(ssp0, period_number*(1/f))
#    eigenValues, eigenVectors = np.linalg.eig(jacobian)
    end = time.time()
#    np.save('/Users/gducci/UCL/Simulations/RefinementTestCase26_100pt/Results/eigen_4states_mod', eigenValues)
    print(end - start)  