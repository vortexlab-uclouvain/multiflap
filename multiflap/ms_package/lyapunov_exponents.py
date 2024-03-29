import numpy as np
from .rk_integrator import rk2, rk3, rk4
import time
import collections
from scipy.integrate import odeint
class LyapunovExponents:

    def __init__(self, x0, n, t_f, ms_object = None):


        """
        Lyapunov exponent class

        Arguments:
            x0: Initial condition for time integration
            n: Number of intermediate points to rescale the trajectory separation
            t_f: Final time of integration
            ms_object: file containind the set of ODEs defined in odes/ folder

        """

        self.n = n
        self.t_f= t_f
        self.ms_object = ms_object
        self.x0 = x0
        self.delta_time = self.t_f/self.n
        #self.d = self.x0 * 1e-5
    def get_lyapunov_exponent(self):

        """ Calculating the largest Lyapunov exponent numerically

                Implementation of the leading Lyapunov exponent obtained rescaling
                the flow separation with the initial perturbation in order to keep.

                Output:

                    lambda_t: Leading Lyapunov exponent

        """
        lambda_t = 0

        # remove x0 transient and let it fall in the attraction domain
        time_array = np.linspace(0, 200, 1000)
        x = odeint(self.ms_object.dynamics, self.x0, time_array)
        x = x[-1,:]
        # perturbation of the initial value
        d = x*1e-9
        x_pert = x + d
        norm_d = np.linalg.norm(d)
        history = np.zeros(self.n)

        lambda_local = np.zeros(self.n)
        for i in range (self.n):
            t0 =i*self.delta_time
            time_array = np.linspace(t0, t0 + self.delta_time, 5000)
            integration_unp = odeint(self.ms_object.dynamics, x, time_array)
            integration_pert = odeint(self.ms_object.dynamics, x_pert, time_array)
            fx = integration_unp[-1, :]
            fx_pert = integration_pert[-1, :]
            d_j = fx_pert - fx
            norm_d_j = np.linalg.norm(d_j)
            lambda_local[i] = np.log(norm_d_j/norm_d) 
            x = fx
            x_pert = x + (d_j*(norm_d/norm_d_j))

        lambda_t = np.sum(lambda_local)/(self.n*self.delta_time)
        return lambda_t
