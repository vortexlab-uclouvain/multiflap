import numpy as np
import bird_model as bm
import rk_integrator as rk
import time

class MultipleShooting:

    def __init__(self, x0, M = 2, model = None):
        self.point_number = M
        self.dimension = 4  # number of states of the ODE system
        self.x0 = x0        # first initial guess
        self.model = model
        self.period = 1/(model.frequency)
        self.time_steps = int(100/(self.point_number - 1))
        self.tau = (self.period)/(self.point_number-1)

    def get_mappedpoint(self,x0, t0, deltat):

        t_final = t0 + deltat     # Final time

        time_array = np.linspace(t0, t_final, self.time_steps)  # Time array for solution
        solution = rk.rk2(self.model.dynamics, x0, time_array)
    #    sspSolution = ode.solve_ivp(birdEqn_py, [tInitial, tFinal], ssp0,'LSODA', max_step = deltat/Nt)
    #    sspSolution = (sspSolution.y).T
        mapped_point = solution[-1, :]  # Read the final point to sspdeltat
        return mapped_point, solution

    def get_initial_guess(self):

        guessed_points = np.zeros([self.point_number,self.dimension])

        # Type here the first guessed point, that will be used to calculate the other (point_number-1) points

        guessed_points[0,0:] = self.x0

        # Automatic routine to extract the remaining point_number-1 points for the flow.
        # Note this is not always the best way to guess points
        for i in range (1,self.point_number):
            [guessed_points[i,0:], _] = self.get_mappedpoint(guessed_points[i-1,:], (i-1)*self.tau, self.tau)
        return guessed_points


    def jacobian_ode(self, x0_jacobian, t):
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
        x0 = x0_jacobian[0:self.dimension]  # First three (or d in general) elements form the original state
                                # space vector
        J = x0_jacobian[self.dimension:].reshape((self.dimension, self.dimension))  # Last nine elements corresponds to
                                             # the elements of Jacobian.
        #We used numpy.reshape function to reshape d^2 self.dimensionensional vector which
        #hold the elements of Jacobian into a dxd matrix.
        #See
        #http://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
        #for the reference for numpy.reshape function

        velJ = np.zeros(np.size(x0_jacobian))  # Initiate the velocity vector as a
                                               # vector of same size with
                                               # x0_jacobian (d + d^2)
    #    velJ[0:self.dimension] = birdEqn_py(t, ssp)
        velJ[0:self.dimension] = self.model.dynamics(x0, t)


        #Last dxd elements of the velJ are determined by the action of
        #stability matrix on the current value of the Jacobian:

        velTangent = np.dot(self.model.get_stability_matrix(x0, t), J)  # Velocity matrix (calculated in the ssp
                                                      # points)  for the tangent space

        velJ[self.dimension:] = np.reshape(velTangent, self.dimension**2)  # Another use of numpy.reshape, here
                                              # to convert from dxd to d^2

        return velJ

    def get_jacobian_analytical(self, x0, initial_time, integration_time):

        """
        Jacobian function for the trajectory started on ssp, evolved for time t

        Inputs:
        ssp: Initial state space point. dx1 NumPy array: ssp = [x, y, z]
        t: Integration time
        Outputs:
        J: Jacobian of trajectory f^t(ssp). dxd NumPy array
        """
        #Hint: See the Jacobian calculation in CycleStability.py
        Jacobian0 = np.identity(self.dimension) #Initial condition for Jacobian matrix
    #    Jacobian0[0,1] = 0.1
        sspJacobian0 = np.zeros(self.dimension + self.dimension ** 2)  # Initiate
        sspJacobian0[0:self.dimension] = x0  # First 3 elemenets
        sspJacobian0[self.dimension:] = np.reshape(Jacobian0, self.dimension**2)  # Remaining 9 elements
        #print(sspJacobian0)
        t_Final = initial_time + integration_time  # Final time
        Nt = 25  # Number of time points to be used in the integration
        tArray = np.linspace(initial_time, t_Final, Nt)  # Time array for solution
        start_jac = time.time()
    #    sspJacobianSolution = ode.solve_ivp(JacobianVelocity,[t_initial, t_Final], sspJacobian0, 'RK45')
        sspJacobianSolution = rk.rk2(self.jacobian_ode, sspJacobian0, tArray)
        end_jac = time.time()
        print("Jacobian time ", (end_jac-start_jac))

    #    sspJacobianSolution = sspJacobianSolution.y.T
        #Read the Jacobian for the periodic orbit:
        J = sspJacobianSolution[-1, self.dimension:].reshape((self.dimension, self.dimension))
        return J


    def get_jacobian_numerical(self, x0, initial_time, integration_time):

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

        jacobian = np.zeros((self.dimension,self.dimension))

        # -------------------------------------------------------------------------
        # Set the numerical perturbation over the direction of the flow
        # -------------------------------------------------------------------------

        epsilon = 1e-3


        # -------------------------------------------------------------------------
        # Finite difference scheme for jacobian evaluation
        # -------------------------------------------------------------------------

        print("... Running jacobian Function")
        for j in range (self.dimension):
            perturbation = np.zeros(self.dimension)

            perturbation[j] = perturbation[j] + x0[j]*epsilon
            x0_pert = x0 + perturbation
            [vel, _] = self.get_mappedpoint(x0, initial_time, integration_time)
            [vel_pert, _] =  self.get_mappedpoint(x0_pert, initial_time, integration_time)
            for i in range (self.dimension):
                jacobian[i,j] = (vel_pert[i] - vel[i])/perturbation[j]


        print("... jacobian Calculated")
        return jacobian

    def get_df_vector(self, x0):
        """


        Similar to shootingMatrix(). Only form the difference matrix
        Parameter: the same as shootingMatrix()
        return: the different vector
        """
        # next line could be omitted and by calling self.dimension and self.point_number.
        # implement change
        M, N = x0.shape
        #xx = states_stack
        dF = np.zeros(M*N)
        #self.tau = nstp*dt
        x_m = x0[-1,:]
        x_0 = x0[0,:]
        for j in range(0, M - 1):
            x_start = x0[j,:]
            fx_start, complete_solution = self.get_mappedpoint(x_start, j*self.tau, self.tau)
            x_end = x0[j+1,:]
            dF[N*j: N*(j+1)] = -(fx_start - x_end)

        dF[-N:] = -(x_m - x_0)

        return dF

    def get_ms_scheme(self, x0):


        """
        Initialization of DF matrix and dF vector
        """
        # The dimension of the MultiShooting matrix is (NxM,NxM)
        DF = np.zeros([self.dimension*self.point_number, self.dimension*(self.point_number)])

        # The dimension of the error vector is (self.dimensionxself.point_number)
        dF = np.zeros(self.dimension*self.point_number)


        # Last guessed point x_m
        x_m = x0[-1,:]

        # starting guessed point x_0
        x_0 = x0[0,:]

        # Routine to fill the rest of the scheme
        complete_solution = []
        for i in range(0, self.point_number - 1):
            x_start = x0[i,:]
            x_end = x0[i+1,:]
            jacobian = self.get_jacobian_analytical(x_start, i*self.tau, self.tau)
            fx_start, trajectory_points = self.get_mappedpoint(x_start, i*self.tau, self.tau)
            complete_solution.append(trajectory_points)
            DF[(i*self.dimension):self.dimension+(i*self.dimension), (i*self.dimension)+self.dimension:2*self.dimension+(i*self.dimension)] = -np.eye(self.dimension)
            DF[(i*self.dimension):self.dimension+(i*self.dimension), (i*self.dimension):(i*self.dimension)+self.dimension] = jacobian
            dF[(i*self.dimension):self.dimension+(i*self.dimension)] = -(fx_start - x_end)

        trajectory = np.asanyarray(complete_solution)
        # Last block of the scheme
        DF[-self.dimension:, 0:self.dimension] = -np.eye(self.dimension)
        DF[-self.dimension:, -self.dimension:] = np.eye(self.dimension)
        dF[-self.dimension:] = -(x_m - x_0)
    #    print(dF)
        print("Error vector is", dF)
        return DF, dF, trajectory
