import numpy as np
import bird_model as bm
import rk_integrator as rk

class MultipleShooting:

    def __init__(self, x0, model = None):
        self.point_number = 2
        self.dimension = 4  # number of states of the ODE system
        self.x0 = x0        # first initial guess
        self.model = model
        self.period = 1/(model.frequency)
        self.time_steps = 100/(self.point_number - 1)
    def get_mappedpoint(self,x0, t0, deltat, time_steps):

        t_final = t0 + deltat     # Final time

        time_array = np.linspace(t0, t_final, time_steps)  # Time array for solution
        solution = rk.rk2(self.model.dynamics, x0, time_array)
    #    sspSolution = ode.solve_ivp(birdEqn_py, [tInitial, tFinal], ssp0,'LSODA', max_step = deltat/Nt)
    #    sspSolution = (sspSolution.y).T
        mapped_point = solution[-1, :]  # Read the final point to sspdeltat
        return mapped_point, solution

    def get_initial_guess(self):

        guessed_points = np.zeros([self.point_number,self.dimension])
        tau = (self.period)/(self.point_number-1)

        # Type here the first guessed point, that will be used to calculate the other (point_number-1) points

        guessed_points[0,0:] = self.x0

        # Automatic routine to extract the remaining point_number-1 points for the flow.
        # Note this is not always the best way to guess points
        for i in range (1,self.point_number):
            [guessed_points[i,0:], _] = self.get_mappedpoint(guessed_points[i-1,:], (i-1)*tau, tau, 30)
        return guessed_points

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
            [vel, _] = self.get_mappedpoint(x0, initial_time, integration_time, time_steps)
            [vel_pert, _] =  self.get_mappedpoint(x0_pert, initial_time, integration_time, time_steps)
            for i in range (self.dimension):
                jacobian[i,j] = (vel_pert[i] - vel[i])/perturbation[j]


        print("... jacobian Calculated")
        return jacobian

    def get_ms_scheme(self, x0, tau, **optional):

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
            Jacob = self.get_jacobian_numerical(x_start, i*tau, tau, **optional)
            fx_start, trajectory_points = self.get_mappedpoint(x_start, i*tau, tau, self.time_steps, **optional)
            complete_solution.append(trajectory_points)
            DF[(i*self.dimension):self.dimension+(i*self.dimension), (i*self.dimension)+self.dimension:2*self.dimension+(i*self.dimension)] = -np.eye(self.dimension)
            DF[(i*self.dimension):self.dimension+(i*self.dimension), (i*self.dimension):(i*self.dimension)+self.dimension] = Jacob
            dF[(i*self.dimension):self.dimension+(i*self.dimension)] = -(fx_start - x_end)

        trajectory = np.asanyarray(complete_solution)
        # Last block of the scheme
        DF[-self.dimension:, 0:self.dimension] = -np.eye(self.dimension)
        DF[-self.dimension:, -self.dimension:] = np.eye(self.dimension)
        dF[-self.dimension:] = -(x_m - x_0)
    #    print(dF)
        print("Error vector is", dF)
        return DF, dF, trajectory
