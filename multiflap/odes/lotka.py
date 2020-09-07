import numpy as np
from scipy.integrate import odeint
class lotka:
    def __init__(self, a=None, b=None, c=None, d=None, e=None, q=None, p=None):

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.q = q
        self.p = p
        self.dimension = 2
    def dynamics(self, x0, t):

        """ODE system
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values
            t: time

        Outputs:
            x_dot: velocity vector
        """
        x, y= x0
        dxdt = 1.5*x - x*y
        dydt = x*y - 3*y

        vel_array = np.array([dxdt, dydt])
        return vel_array
