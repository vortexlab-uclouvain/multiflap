import numpy as np  # Import NumPy
from scipy.integrate import odeint  # Import odeint
from RungeKutta import RK4 
dim = 4
g = 9.81

def Velocity(ssp, t):
    """
    Velocity function for the Rossler flow

    Inputs:
    ssp: State space vector. dx1 NumPy array: ssp=[x, y, z]
    t: Time. Has no effect on the function, we have it as an input so that our
       ODE would be compatible for use with generic integrators from
       scipy.integrate

    Outputs:
    vel: velocity at ssp. dx1 NumPy array: vel = [dx/dt, dy/dt, dz/dt]
    """
    u, w, theta, q = ssp  # Read state space points
    # Rossler flow equations:
    dudt = -q*w - g*np.sin(theta) 
    dwdt = q*u + g*np.cos(theta)
    dqdt = 0
    dthetadt = q
    # Collect Equations of motion in a single NumPy array:
    vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector

    return vel

#def StabilityMatrix(ssp):
#    """
#    Stability matrix for the Rossler flow
#
#    Inputs:
#    ssp: State space vector. dx1 NumPy array: ssp = [x, y, z]
#    Outputs:
#    A: Stability matrix evaluated at ssp. dxd NumPy array
#       A[i, j] = del Velocity[i] / del ssp[j]
#    """
#
#    x, y, z = ssp  # Read state space points
#
#    A = np.array([[0, -1, -1],
#                  [1, a, 0],
#                  [z, 0, x-c]], float)  # COMPLETE THIS LINE
#    return A

def Flow(ssp0, deltat):
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
    tInitial = 0  # Initial time
    tFinal = deltat  # Final time
    Nt = 1000 # Number of time points to be used in the integration

    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array for solution

    sspSolution = RK4(Velocity, ssp0, tArray)
    sspdeltat = sspSolution[-1, :]  # Read the final point to sspdeltat
    return sspdeltat


#def JacobianVelocity(sspJacobian, t):
#    """
#    Velocity function for the Jacobian integration
#
#    Inputs:
#    sspJacobian: (d+d^2)x1 dimensional state space vector including both the
#                 state space itself and the tangent space
#    t: Time. Has no effect on the function, we have it as an input so that our
#       ODE would be compatible for use with generic integrators from
#       scipy.integrate
#
#    Outputs:
#    velJ = (d+d^2)x1 dimensional velocity vector
#    """
#    
#    # Prendo il vettore sspJacobian (d + d^2) e la parte d la riempio, mentre 
#    # la parte d^2 la trasformo nel Jacobiano (d x d)
#    ssp = sspJacobian[0:dim]  # First three (or d in general) elements form the original state
#                            # space vector
#    J = sspJacobian[dim:].reshape((dim, dim))  # Last nine elements corresponds to
#                                         # the elements of Jacobian.
#    #We used numpy.reshape function to reshape d^2 dimensional vector which
#    #hold the elements of Jacobian into a dxd matrix.
#    #See
#    #http://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
#    #for the reference for numpy.reshape function
#
#    velJ = np.zeros(np.size(sspJacobian))  # Initiate the velocity vector as a
#                                           # vector of same size with
#                                           # sspJacobian (d + d^2)
#    velJ[0:dim] = Velocity(ssp, t)
#    
#   
#    #Last dxd elements of the velJ are determined by the action of
#    #stability matrix on the current value of the Jacobian:
#
#    velTangent = np.dot(StabilityMatrix(ssp), J)  # Velocity matrix (calculated in the ssp
#                                                  # points)  for the tangent space
#
#    velJ[dim:] = np.reshape(velTangent, dim**2)  # Another use of numpy.reshape, here
#                                          # to convert from dxd to d^2
#
#    return velJ


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
    sspJacobian0 = np.zeros(dim + dim ** 2)  # Initiate
    sspJacobian0[0:dim] = ssp  # First 3 elemenets
    sspJacobian0[dim:] = np.reshape(Jacobian0, dim**2)  # Remaining 9 elements
    #print(sspJacobian0)
    tInitial = 0  # Initial time
    tFinal = t  # Final time
    Nt = 100000  # Number of time points to be used in the integration
    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array for solution

    sspJacobianSolution = odeint(JacobianVelocity, sspJacobian0, tArray)
    x1t = sspJacobianSolution[:, 0]  # Read x(t)
    x2t = sspJacobianSolution[:, 1]  # Read y(t)
    x3t = sspJacobianSolution[:, 2]  # Read z(t)
    x4t = sspJacobianSolution[:, 3]  # Read z(t)

    #Read the Jacobian for the periodic orbit:
    J = sspJacobianSolution[-1, dim:].reshape((dim, dim))
    return J


# This block will be evaluated if this script is called as the main routine
# and will be ignored if this file is imported from another script.
# To be used as a debug/test.
    
if __name__ == "__main__":
   ssp = [0.00702625,	-0.035131,	0.035131]
   t = 10*6.881088454743404
   tArray = np.linspace(0,t,10000)
   test = RK4(Velocity, ssp,tArray)
   print(test)
