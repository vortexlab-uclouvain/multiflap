import numpy as np


class Jaschinski:
    def __init__(self, V = 20.*math.sqrt(5),
                  c_x = 4e3*12.5,
                  c_y = 300*12.5,
                  c_z = 6.25e5,
                  m = 1005.,
                  In = 572.,
                  In_y = 95,
                  a_0 = 0.753,
                  r_0 = 0.5,
                  delta_0 = 0.0262,
                  R_R = 0.3,
                  g = 9.81
                  L = 0.95,
                  G = 7.92e10,
                  sigma = 0.28,
                  mu = 0.12
                  y0 = 0.05):

        self.V = V
        self.c_x = c_x
        self.c_y = c_y
        self.c_z = c_z
        self.m = m
        self.In = In
        self.In_y = In_y
        self.a_0 = a_0
        self.r_0 = r_0
        self.delta_0 = delta_0
        self.R_R = R_R
        self.g = g
        self.L = L
        self.G = G
        self.sigma = sigma
        self.mu = mu
        self.dimension = 4
        self.Gamma = self.delta_0/(self.a_0 - self.r_0*sef/delta_0)
        self.khi = self.Gamma*self.a_0/self.delta_0
        self.c_0 = self.delta_0**2 / self.Gamma
        self.b_0 = 2*self.Gamma + self.Gamma**2*(self.R_R+self.r_0)
        self.y0 = y0
        self.sign_init = 1
        self.y = self.y0*self.sign_init
        self.y_d = 0.0* self.sign_init
        self.psi = 0.*self.sign_init
        self.psi_d = 0.0*self.sign_init
        self.N_L = self.m*self.g/2
        self.N_R = self.m*self.g/2
        self.a_L = self.a_0 + self.a_0*self.Gamma/self.delta_0 * self.y + self.R_R*self.Gamma*self.y
        self.a_R = self.a_0 - self.a_0*self.Gamma/self.delta_0 * self.y - self.R_R*self.Gamma*self.y
        self.r_L = self.r_0 - self.a_0*self.Gamma*self.y - self.R_R*self.delta_0*self.Gamma*self.y
        self.r_R = self.r_0 - self.a_0*self.Gamma*self.y + self.R_R*self.delta_0*self.Gamma*self.y
        self.xi_L = -self.r_L*self.psi*delta_0
        self.xi_R = self.r_R*self.psi*delta_0

    def HertzKalker(wheel_radius, rail_radius):
        A = 1/wheel_radius;
        B = 1/rail_radius;

        # theta
        theta = math.acos(abs((B-A)/(B+A)))

        if(theta < 0.0 or theta>= math.pi/2):{
            print("\t >Hertz> theta ({:.3f}) is outside validity range of m,n.\n".format(theta))
        }

        m_ratio = 0.84 + 0.3981e-3 * (4.05 - theta)** 6.709;
        n_ratio = 0.22 + 0.1483    * (0.88 + theta)** 1.845;

        ratio = n_ratio/m_ratio;

        C11 = 3.2893 + 0.9750 / ratio - 0.0120 / (ratio **2)
        C22 = 2.4014 + 1.3179 / ratio - 0.0200 / (ratio **2)
        C23 = 0.4147 + 1.0184 / ratio + 0.0565 / (ratio **2) - 0.0013 / (ratio **3)

        return C11, C22, C23

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
        # Position of contact points (Eq 165 p 111)
        a_L = self.a_0 + self.a_0*self.Gamma/self.delta_0 * self.y + self.R_R*self.Gamma*self.y
        a_R = a_0 - a_0*Gamma/delta_0 * y - R_R*Gamma*y
        r_L = r_0 - a_0*Gamma*y - R_R*delta_0*Gamma*y
        r_R = r_0 + a_0*Gamma*y + R_R*delta_0*Gamma*y
        xi_L = -r_L*psi*delta_0
        xi_R =  r_R*psi*delta_0

        # other variables
        # Eq 57 p30
        s = khi*r_0/2./a_0 - delta_0*khi**2*In/(2.*a_0**2*m) + delta_0/2.*khi
        # alpha defined just after s (p30)
        alpha = 0#delta_0*psi

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



        # Average radius of contact ellipse (p33, Eq.213-p139)
        c_L = 19.027e-5*(N_L * 2.*r_0*R_R/(r_0+R_R))**(1./3.)
        c_R = 19.027e-5*(N_R * 2.*r_0*R_R/(r_0+R_R))**(1./3.)

        # Kalker coef
        C11_L, C22_L, C23_L = HertzKalker(r_L, R_R)
        C11_R, C22_R, C23_R = HertzKalker(r_R, R_R) 


        # Abreviation of  Kalker coef (p30)
        f11_L = G*c_L**2*C11_L
        f11_R = G*c_R**2*C11_R
        f22_L = G*c_L**2*C22_L
        f22_R = G*c_R**2*C22_R
        f23_L = G*c_L**3*C23_L
        f23_R = G*c_R**3*C23_R

        # Kalker forces (Eq 63 p 33)
        T_x_L = -mu*N_L * math.tanh(f11_L/(mu*N_L)*( a_0*Gamma/r_0 *y + a_0/V*psi_d))
        T_x_R = -mu*N_R * math.tanh(f11_R/(mu*N_R)*(-a_0*Gamma/r_0 *y - a_0/V*psi_d))
        T_y_L = -mu*N_L * math.tanh(f22_L/(mu*N_L)*(-psi+khi/V*y_d) + f23_L/(mu*N_L)*( delta_0/r_0 + psi_d/V))
        T_y_R = -mu*N_R * math.tanh(f22_R/(mu*N_R)*(-psi+khi/V*y_d) + f23_R/(mu*N_R)*(-delta_0/r_0 + psi_d/V))

        T_L = math.sqrt(T_x_L**2+T_y_L**2)
        T_prime_L = mu*N_L
        if T_L>T_prime_L:
            T_x_L = T_x_L*T_prime_L/T_L
            T_y_L = T_y_L*T_prime_L/T_L

        T_R = math.sqrt(T_x_R**2+T_y_R**2)
        T_prime_R = mu*N_R
        if T_R>T_prime_R:
            T_x_R = T_x_R*T_prime_R/T_R
            T_y_R = T_y_R*T_prime_R/T_R

        # Eq 72 (p38) in form M*x = C + Fk1 + Kk2
        C = np.zeros(4)
        C[0] =  In_y * Gamma*V/khi/r_0 * psi_d - m*g* b_0/khi * y - 2*c_y/khi * y
        C[1] = -In_y * Gamma*V/r_0     * y_d   + m*g* psi* c_0 - 2*c_x* L**2 *psi
        C[2] =  In_y * khi*V/(2.*r_0*a_0) * psi_d + m*g/2. - khi**2 * m*g/2./a_0 *y + c_y *delta_0*khi**2*In/(a_0**2*m)*y - c_z *m* L**2 *delta_0/2/In *y
        C[3] = -In_y * khi*V/(2.*r_0*a_0) * psi_d + m*g/2. + khi**2 * m*g/2./a_0 *y - c_y *delta_0*khi**2*In/(a_0**2*m)*y + c_z *m* L**2 *delta_0/2/In *y

        Fk1 = np.zeros(4)
        Fk1[0] = T_y_L + T_y_R + psi * (T_x_L + T_x_R)
        Fk1[1] = a_0*(T_x_L - T_x_R) + khi*y*(T_x_L + T_x_R)
        Fk1[2] =  s*(T_y_L + T_y_R ) + delta_0/2.*(T_y_L - T_y_R )
        Fk1[3] = -s*(T_y_L + T_y_R ) + delta_0/2.*(T_y_L - T_y_R )

        Fk2 = np.zeros(4)
        Fk2[0] = 0.
        Fk2[1] = -alpha*r_0*(T_y_L - T_y_R )
        Fk2[2] = alpha * khi*(T_x_L + T_x_R)*(.5 - khi*In/(2*m*a_0**2)) + alpha*(T_x_L - T_x_R)*(.5 - m*a_0**2/2./In)
        Fk2[3] = alpha * khi*(T_x_L + T_x_R)*(.5 - khi*In/(2*m*a_0**2)) + alpha*(T_x_L - T_x_R)*(.5 - m*a_0**2/2./In)

        # compute right and side
        rhs = C + Fk1 + Fk2

        y_dd = rhs[0]/(m/khi)
        psi_dd = rhs[1]/In
        N_L = rhs[2]
        N_R = rhs[3]

        # print("N_L={} N_R={}".format(N_L, N_R))
        
        # print(C)
        # print(Fk1)
        # print(Fk2)
        # print("T_x_L={} T_x_R={} T_y_L={}   T_y_R={}".format(T_x_L, T_x_R, T_y_L, T_y_R))
        # print("N_L={} N_R={} y_dd={}   psi_dd={}".format(N_L, N_R, y_dd, psi_dd))
        # print(" ")
        return y_dd, psi_dd
