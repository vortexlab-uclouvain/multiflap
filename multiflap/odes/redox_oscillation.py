import numpy as np

"""
Example case adopted from:

    A Robust Model for Circadian Redox Oscillations,
    del Olmo, M.; Kramer, A.; Herzel, H., Int. J. Mol. Sci., 2019, 20, 2368
    https://www.mdpi.com/1422-0067/20/9/2368

"""


class RedoxModel:
    def __init__(self, a=1000, b=2, c=10000, d=0.2, e=0.1, q=0.1, p=1):

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.q = q
        self.p = p
        self.dimension = 4
    def dynamics(self, x0, t):

        """ODE system
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values
            t: time

        Outputs:
            x_dot: velocity vector
        """
        D1, D2, R, A = x0
        dD1_dt = self.p - self.a*A*D1 - self.d*D1
        dD2_dt = self.d*D1 - self.e*D2
        dR_dt = self.e*D2 - self.q*R
        dA_dt = self.b*(1-A)*R - self.a*A*D1

        vel_array = np.array([dD1_dt, dD2_dt, dR_dt, dA_dt], float)
        return vel_array



    def get_stability_matrix(self, x0, t):

        """
        Stability matrix of the ODE system for the longitudinal plane

        Inputs:
            x0: initial condition
        Outputs:
            A: Stability matrix evaluated at x0. (dxd) dimension
            A[i, j] = dv[i]/dx[j]
        """
        D1, D2, R, A = x0
        A_matrix = np.array([[-self.d - self.a*A,  0., 0.,-self.a*D1],
                      [self.d,  -self.e, 0., 0.],
                      [0., self.e, -self.q, 0.],
                      [-self.a*A, 0., self.b*(1-A), -self.b*R -self.a*D1]], float)

        return A_matrix
