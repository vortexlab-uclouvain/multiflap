import numpy as np

"""
Case adopted from:

    Optimized shooting method for finding periodic orbits of nonlinear dynamical systems
    Dednam, W and Botha, Andre E
    https://arxiv.org/abs/1405.5347

"""

class Rossler:
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
            x0: initial values
            t: time

        Outputs:
            x_dot: velocity vector
        """
        x, y, z = x0
        dxdt = - y - z
        dydt = x + self.a*y
        dzdt = self.b + z*(x - self.c)
        vel_array = np.array([dxdt, dydt, dzdt], float)
        return vel_array



    def get_stability_matrix(self, x0, t):

        """
        Stability matrix of the ODE system

        Inputs:
            x0: initial condition 
        Outputs:
            A: Stability matrix evaluated at x0. (dxd) dimension
            A[i, j] = dv[i]/dx[j]
        """
        x, y, z = x0
        A_matrix = np.array([[0, -1, -1],
                            [1, self.a, 0],
                            [z, 0, x-self.c]], float)

        return A_matrix
