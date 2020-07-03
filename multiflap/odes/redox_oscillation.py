import numpy as np


class RedoxModel:
    def __init__(self, a=1000, b=2, c=10000, d=0.2, e=0.1, q=0.1, p=1):

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.q = q
        self.p = p

    def dynamics(self, x0, t):

        D1, D2, R, A = x0
        dD1_dt = self.p - self.a*A*D1 - self.d*D1
        dD2_dt = self.d*D1 - self.e*D2
        dR_dt = self.e*D2 - self.q*R
        dA_dt = self.b*(1-A)*R - self.a*A*D1

        vel_array = np.array([dD1_dt, dD2_dt, dR_dt, dA_dt], float)
        return vel_array
