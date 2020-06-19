import numpy as np
import FlowFunctions_V1 as func
#from scipy.sparse.linalg import gmres
from numpy.linalg import norm
from numpy import inf
import os
import time

"""
Parameters of the shooting scheme:
    dim = dimension of the ODE system, this is read from the FlowFunction file
    M = The number of guessed points is set here
    N = equivaent of "dim". It is just used for a better readibility of the code
    The dimension of the MultipleShooting matrix is (M x N)
"""

# Dimension of the problem read from FlowFunction

dim = func.dim

# Number of guessed points (integer)
M = 2

# Dimension of the problem under another name
N = dim


time_steps = 100
"""
MultiShootingScheme function

---------------------------------------------------------------------------
---------------------------------------------------------------------------

    Args:
           states_stack:     Coordinates of the points at each time step
           tau:              Time interval between two consecutives points
           maxIter:          Maximum number of iteration, defined in the main
           tol:              Tolerance before exiting the loop

---------------------------------------------------------------------------
---------------------------------------------------------------------------

    Return:
           x:                Updated vector of points.

    At convengence x belong to limit cycle
---------------------------------------------------------------------------
---------------------------------------------------------------------------

    Newton scheme to update the value of the guessed points at
    every step of the iteration

    (0) DF(step k)*dx(step k) = dF(step k)    <--------------
    (1) Evaluation of the error dF(step k)                   |
    (2) Evaluation of the dx vector                          |
                    dx(step k) = [DF(step k)]^-1*dF(step k)  |
    (3) Update of the guessed points                         |
                    x = x(step k) + dx                       |
    (4) back to step 0 until convergence is obtained --------
                            |
                            |
                            V
                            x(step_k)

"""
def MultiShootingScheme(states_stack, tau, maxIter, tol, result_directory, **optional):

    # Step zero of the scheme, storage of guessed points given as an input
    # by the user
    x = 1*(states_stack)
    lam= 1.0
    #x_new = np.copy(x)
#    data_sim_history = (result_directory+"/Simulation_History")
#    os.mkdir(data_sim_history)

    ptlist = [states_stack]
    error = []
    # Start the iterative scheme
    for i in range(maxIter):

        # For the step "i", evaluation of DF matrix and dF vector, calling the
        # MultiShootingMatrix function
        start_timing = time.time()
        print("... Iteration i = " + str(i)+" ... RUNNING with amplitude ")

        DF, dF, complete_solution = MultiShootingMatrix(x, tau, **optional)

        # Evaluation of the maximum error in the dF error vector
        epsilon = max(np.abs(dF))
        error.append(epsilon)
        end_timing = time.time()
        print("... Iteration " + str(i) + " time: " + str(end_timing - start_timing))
        print("... Iteration number i = " + str(i) + " has error: " + str(norm(dF, inf)))

        # -------------------------------------------------------------------------
        # Save Data every iteration
        # -------------------------------------------------------------------------

#        np.save(data_sim_history+"/Checkpoint_Iteration_"+str(i), x)
#        np.save(data_sim_history+"/Error_History", error)


        # Refinement of the solution, by calculatin the dx vector
        if epsilon < tol:
            print("Iteration terminates at error : " + str(norm(dF, inf)))
            print("Calculating Floquet Multipliers")
#            Jacobian = func.JacobianNumerical(x, 0, (M - 1)*tau)
            ptlist.append(1*x)
            Jacobian_semigroup = np.eye(N)
            # =============================================================================
            #    Evaluation of Jacobian, using Semigroup property
            # =============================================================================
            for i in range(0, M - 1):
                Jacobian_semigroup = (DF[(i*N):N+(i*N), (i*N):(i*N)+N]).dot(Jacobian_semigroup)

            break

        else:
            JJ = np.dot(DF.T,  DF)
            JdF = np.dot(DF.T, dF)
            H = JJ + lam* np.diag(np.diag(JJ))
            delta_x = np.linalg.solve(H, JdF)
            #delta_x = np.linalg.solve((DF), dF) # Matrice dx che trovo invertendo il sistema

            # Update of the solution
            for k in range(M):
                x[k,:] = x[k,:] + delta_x[k*N:(k+1)*N]
            dF_new = dFvector(x, tau, **optional)
            if norm(dF_new, inf) < norm(dF, inf):
                #x = x_new;
                tau = tau
                lam= lam / 10.0

            else :
                lam= lam * 10.0
                if lam> 1e10:
                    print ("lam is too large")
            ptlist.append(1*x)

    return x, ptlist, error, complete_solution, Jacobian_semigroup


def dFvector(states_stack, tau, **optional):
    """


    Similar to shootingMatrix(). Only form the difference matrix
    Parameter: the same as shootingMatrix()
    return: the different vector
    """
    M, N = states_stack.shape
    #xx = states_stack
    dF = np.zeros(M*N)
    #tau = nstp*dt
    x_m = states_stack[-1,:]
    x_0 = states_stack[0,:]
    for j in range(0, M - 1):
        x_start = states_stack[j,:]
        fx_start, complete_solution = func.Flow(x_start, j*tau, tau, time_steps, **optional)
        x_end = states_stack[j+1,:]
        dF[N*j: N*(j+1)] = -(fx_start - x_end)

    dF[-N:] = -(x_m - x_0)

    return dF


"""
Function MultiShootingMatrix

    Args:  states_stack point (NxM): contains the points of the orbit at the step k
           dt is frozen from the main
           nstp is frozen from the main

    return: DF matrix, dF vector

    DF * dx = dF

Here is the scheme
---------------------------------------------------------------------------
---------------------------------------------------------------------------

    Multi Shooting Matrix:

        dim(DF) = ((NxM) x (NxM))

                 [  J(x_0, tau)  -I                                  ]
                 [              J(x_1, tau)   -I                     ]
        DF   =   [                .                                  ]
                 [                 .                                 ]
                 [                  .                                ]
                 [                        J(x_{M-1}, tau)       I    ]
                 [ -I ........................................  I    ]

---------------------------------------------------------------------------
---------------------------------------------------------------------------

    Unknown Vector:

        dim(dx) = (NxM)

                [dx_0]
        dx   =  [dx_1]
                [... ]
                [dx_M]

---------------------------------------------------------------------------
---------------------------------------------------------------------------

    Error vector:

        dim(dF) = (NxM)

                [ f(x_0, tau) - x_1     ]
        dF  = - [ f(x_1, tau) - x_2     ]
                [           ...         ]
                [       x_M - x_0       ]

---------------------------------------------------------------------------
---------------------------------------------------------------------------
"""

def MultiShootingMatrix(states_stack, tau, **optional):

    # Time interval between two consecutives points
    #tau = nstp*dt

    """
    Initialization of DF matrix and dF vector
    """
    # The dimension of the MultiShooting matrix is (NxM,NxM)
    DF = np.zeros([N*M, N*(M)])

    # The dimension of the error vector is (NxM)
    dF = np.zeros(N*M)


    # Last guessed point x_m
    x_m = states_stack[-1,:]

    # starting guessed point x_0
    x_0 = states_stack[0,:]

    # Routine to fill the rest of the scheme
    complete_solution = []
    for i in range(0, M - 1):
        x_start = states_stack[i,:]
        x_end = states_stack[i+1,:]
        Jacob = func.Jacobian(x_start, i*tau, tau, **optional)
        fx_start, trajectory_points = func.Flow(x_start, i*tau, tau, time_steps, **optional)
        complete_solution.append(trajectory_points)
        DF[(i*N):N+(i*N), (i*N)+N:2*N+(i*N)] = -np.eye(N)
        DF[(i*N):N+(i*N), (i*N):(i*N)+N] = Jacob
        dF[(i*N):N+(i*N)] = -(fx_start - x_end)

    trajectory = np.asanyarray(complete_solution)
    # Last block of the scheme
    DF[-N:, 0:N] = -np.eye(N)
    DF[-N:, -N:] = np.eye(N)
    dF[-N:] = -(x_m - x_0)
#    print(dF)
    print("Error vector is", dF)
    return DF, dF, trajectory
