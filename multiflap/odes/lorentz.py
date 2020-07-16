import numpy as np


class Lorentz:
    def __init__(self, a=None, b=None, c=None, d=None, e=None, q=None, p=None):

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.q = q
        self.p = p
        self.dimension = 3
    def dynamics(self, x0, t):

        """ODE system
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values [u, w, q, theta]
            t: time

        Outputs:
            x_dot: velocity vector
        """
        x, y, z = x0
        dxdt = self.a*(y-x)
        dydt = x*(self.b-z)-y
        dzdt = x*y - self.c*z 
        vel_array = np.array([dxdt, dydt, dzdt], float)
        return vel_array



    def get_stability_matrix(self, x0, t):

        """
        Stability matrix of the ODE system for the longitudinal plane

        Inputs:
            x0: initial condition [u_0, w_0, q_0, theta_0]
        Outputs:
            A: Stability matrix evaluated at x0. (dxd) dimension
            A[i, j] = dv[i]/dx[j]
        """
        x, y, z = x0
        A_matrix = np.array([[-self.a,  self.a, 0],
                            [(self.b-z), -1, -x],
                            [y, x, -self.c]], float)

        return A_matrix
