import numpy as np
from .rk_integrator import rk2, rk3, rk4
import time
import collections

class MultipleShooting:

    def __init__(self, x0, M = 2, period = 0.25,
                 t_steps = 70, model = None, option_jacobian = 'analytical'):
        self.point_number = M
        self.dim = model.dimension  # number of states of the ODE system
        self.x0 = x0        # first initial guess
#        if not isinstance(model, Model):
 #           raise TypeError("Model missing!")
        self.model = model
        self.period = period
        self.option_jacobian = option_jacobian
        self.t_steps = t_steps
        self.tau = (self.period)/(self.point_number-1)

    def get_mappedpoint(self,x0, t0, deltat):
        """
        Returns the last element of the time integration. It outputs where a
        point x0(t) is mapped after a time deltat.

        Inputs:
            x0: initial value
            t0: initial time (required as the system is non-autonomous)
            deltat: integration_time

        Outputs:
            mapped_point: last element of the time integration
            solution: complete trajectory traced out from x0(t0) for t = deltat


        """
        t_final = t0 + deltat     # Final time

        time_array = np.linspace(t0, t_final, self.t_steps)
        rk_solution = rk4(self.model.dynamics, x0, time_array)
    #    sspSolution = ode.solve_ivp(birdEqn_py, 
                        #[tInitial, tFinal], ssp0,'LSODA', max_step = deltat/Nt)
    #    sspSolution = (sspSolution.y).T
        solution = rk_solution.x
        mapped_point = solution[-1, :]  # Read the final point to sspdeltat
        return mapped_point, rk_solution

    def get_initial_guess(self):

        """
        Starting from x0, it spreads the other M-1 points in the phase space
        just by time integration. The mismatch at the first iteration therefore
        will be just between the point 0 and the point M.

        Input:
            point_number: number of points for the multiple-shooting algorithm
            dimension: dimension of the ODE system
            x0: initial value at time t

        Output:
            guessed_pts: guessed points spread in the phase space


        """
        guessed_pts = np.zeros([self.point_number,self.dim])


        guessed_pts[0,0:] = self.x0

        # Automatic routine to extract the remaining point_number-1 points
        # Note this is just one simple model to spread points. Other methods
        # may be used.
        for i in range (1,self.point_number):
            [guessed_pts[i,0:], _] = self.get_mappedpoint(guessed_pts[i-1,:],
                                                             (i-1)*self.tau,
                                                             self.tau)

        return guessed_pts


    def jacobian_ode(self, x0_jacobian, t):
        """
        Set up the additional ODE system (d + d^2) to evaluate analytically the
        Jacobian matrix.
        This function is used to unpack the Jacobian elements to solve

        \dot{J} = AJ

        It reshapes Equation (7.18) of Seydel's book "Practical Bifurcation and
        Stability Analysis"
        in order get the components of the Jacobian via numerical integration
        Inputs:
            x0_jacobian: (d+d^2) initial values
                         state space itself and the tangent space
            t: time.

        Outputs:
            jac_elements_ODE = (d+d^2) dimensional velocity vector

        """



        x0 = x0_jacobian[0:self.dim]
        J = x0_jacobian[self.dim:].reshape((self.dim, self.dim))

        jac_elements_ODE = np.zeros(np.size(x0_jacobian))

        jac_elements_ODE[0:self.dim] = self.model.dynamics(x0, t)


        #Last dxd elements of the velJ are determined by the action of
        #stability matrix on the current value of the Jacobian:

        velTangent = np.dot(self.model.get_stability_matrix(x0, t), J)

        # shape a back a (dxd) array in a d^2 matrix

        jac_elements_ODE[self.dim:] = np.reshape(velTangent,
                                                       self.dim**2)

        return jac_elements_ODE

    def get_jacobian_analytical(self, x0, initial_time, integration_time):

        """
        Return the Jacobian (or Monodromy Matrix) of the flow, starting from x0
        and integrated for a time "integration_time".

        It solves numerically the (d + d^2) ODE system.

        Reference: Equation (7.18) Seydel's book "Practical Bifurcation and
        Stability Analysis".

        Inputs:
            x0 : Initial point of the phase space. len(x0) = dimension
            initial_time: initial time needed as the system is non-autonomous
            integration_time: integration time
        Outputs:
            J: Jacobian (Monodromy Matrix) of flow from t -> t+integration_time
        """

        # Initial conditions (ic) for Jacobian matrix (see 7.18 Seydel)

        jac_ic = np.identity(self.dim)

        jacODE_ic = np.zeros(self.dim + self.dim ** 2)

        jacODE_ic[0:self.dim] = x0

        jacODE_ic[self.dim:] = np.reshape(jac_ic, self.dim**2)


        t_Final = initial_time + integration_time
        Nt = 30 #50000  # interval discretization for computing the integration

        tArray = np.linspace(initial_time, t_Final, Nt)

        start_jac = time.time()

    #   jac_elements_solution = ode.solve_ivp(jacobian_ode,[t_initial, t_Final],
                                #jacODE_ic, 'RK45')

        rk_jac_elements_solution = rk2(self.jacobian_ode, jacODE_ic, tArray)

        end_jac = time.time()
        print("Jacobian time ", (end_jac-start_jac))

        jac_elements_solution = rk_jac_elements_solution.x
    #    jac_elements_solution = jac_elements_solution.y.T
        # Pack back the jacobian in matrix:
        J = jac_elements_solution[-1, self.dim:].reshape((self.dim,
                                                                self.dim))

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
        # ---------------------------------------------------------------------
        #  Initialization of the Jacobian Matrix
        # ---------------------------------------------------------------------

        jacobian = np.zeros((self.dim,self.dim))

        # ---------------------------------------------------------------------
        # Set the numerical perturbation over the direction of the flow
        # ---------------------------------------------------------------------

        epsilon = 1e-3


        # ---------------------------------------------------------------------
        # Finite difference scheme for jacobian evaluation
        # ---------------------------------------------------------------------

        print("... Running jacobian Function")
        for j in range (self.dim):
            perturbation = np.zeros(self.dim)

            perturbation[j] = perturbation[j] + x0[j]*epsilon

            x0_pert = x0 + perturbation

            [vel, _] = self.get_mappedpoint(x0,
                                            initial_time,
                                            integration_time)

            [vel_pert, _] =  self.get_mappedpoint(x0_pert,
                                                  initial_time,
                                                  integration_time)


            for i in range (self.dim):

                jacobian[i,j] = (vel_pert[i] - vel[i])/perturbation[j]


        print("... jacobian Calculated")
        return jacobian

    def get_error_vector(self, x0):
        """
        Returns the right-hand side of the multiple-shooting scheme
        Error vector -(f(x_i) - x_{i+1})
        dim: n x M

        Inputs:
            x0: array of points of multiple-shooting scheme

        Output:
            error: error vector, rhs of the multiple shooting scheme

        """
        # implement change
        m, n = x0.shape
        #xx = states_stack
        error = np.zeros(m*n)
        #self.tau = nstp*dt
        x_m = x0[-1,:]
        x_0 = x0[0,:]
        LimitCycle = collections.namedtuple('LimitCycle',['space', 'time'])
        full_trajectory = np.zeros([(self.t_steps*(self.point_number-1)), n])
        time = np.zeros(self.t_steps*(self.point_number-1))
        steps = self.t_steps # just to reduce length lines
        for j in range(0, m - 1):
            x_start = x0[j,:]
            fx_start, trajectory_tuple = self.get_mappedpoint(x_start,
                                                               j*self.tau,
                                                               self.tau)
            trajectory = trajectory_tuple.x
            delta_time = trajectory_tuple.t

            full_trajectory[j*(steps):j*(steps)+(steps),:] = trajectory[0:,:]

            time[j*(steps):j*(steps)+(steps)] = delta_time
            x_end = x0[j+1,:]
            error[n*j: n*(j+1)] = -(fx_start - x_end)

        error[-n:] = -(x_m - x_0)
        solution_tuple = LimitCycle(full_trajectory, time)
        return error, solution_tuple

    def get_ms_scheme(self, x0):


        """
        Returns the multiple-shooting scheme to set up the equation:
        MS * DX = E

        Inputs:
            x0: array of m-points of the multiple-shooting scheme

        Outputs (M = points number, N = dimension):
            MS = Multiple-Shooting matrix dim(NxM, NxM)
            error = error vector dim(NxM),
            trajectory_tuple = trajectory related to the initial value and
                                time
        -----------------------------------------------------------------------

        Multiple-Shooting Matrix (MS):

            dim(MS) = ((NxM) x (NxM))

                     [  J(x_0, tau)  -I                                  ]
                     [              J(x_1, tau)   -I                     ]
            MS   =   [                .                                  ]
                     [                 .                                 ]
                     [                  .                                ]
                     [                        J(x_{M-1}, tau)       I    ]
                     [ -I ........................................  I    ]


        Unknown Vector:                 Error vector:

        dim(DX) = (NxM)                 dim(E) = (NxM)

                [DX_0]                          [ f(x_0, tau) - x_1 ]
        DX   =  [DX_1]                  E  =  - [ f(x_1, tau) - x_2 ]
                [... ]                          [       ...         ]
                [DX_M]                          [       x_M - x_0   ]

        """
        # The dimension of the MultiShooting matrix is (NxM,NxM)
        MS = np.zeros([self.dim*self.point_number,
                       self.dim*(self.point_number)])

        # Routine to fill the rest of the scheme
        #complete_solution = []
        for i in range(0, self.point_number - 1):
            x_start = x0[i,:]
            if self.option_jacobian == 'analytical':
                jacobian = self.get_jacobian_analytical(x_start, i*self.tau,
                                                    self.tau)
            if self.option_jacobian == 'numerical':
                jacobian = self.get_jacobian_numerical(x_start, i*self.tau,
                                                        self.tau)


            MS[(i*self.dim):self.dim+(i*self.dim),
               (i*self.dim)+self.dim:2*self.dim+(i*self.dim)]=-np.eye(self.dim)


            MS[(i*self.dim):self.dim+(i*self.dim),
               (i*self.dim):(i*self.dim)+self.dim] = jacobian


        # Last block of the scheme
        MS[-self.dim:, 0:self.dim] = -np.eye(self.dim)
        MS[-self.dim:, -self.dim:] = np.eye(self.dim)

        [error, trajectory_tuple] = self.get_error_vector(x0)
        print("Error vector is", error)
        return MS, error, trajectory_tuple
