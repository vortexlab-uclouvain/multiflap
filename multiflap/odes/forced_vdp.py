import numpy as np

"""
Example case adopted from:

    Practical Bifurcation and Stability Analysis, page 347
    Seydel R.
    Forced Van der Pol oscillator


    gamma = 0.5 ; lambda = 1.8
    p:q = 1:3 (from Seydel results). Meaning the period T=10.4719755

    This example does not use on purpose the stability matrix, in order to
    calculate the Jacoian numerically only.
"""

class ForcedVanDerPol:
    def __init__(self, lam=1.8 , gamma = 0.5, delta=4):
        self.lam = lam
        self.gamma = gamma
        self.delta = delta
        self.dimension=2
    def dynamics(self, x0, t):

        """ODE system
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values
            t: time

        Outputs:
            x_dot: velocity vector
        """
        y1, y2 = x0
        dy1_dt = y2
        dy2_dt = self.delta*(1 - y1**2)*y2 -y1 + self.gamma*(np.cos(self.lam*t))

        vel_array = np.array([dy1_dt, dy2_dt], float)
        return vel_array
