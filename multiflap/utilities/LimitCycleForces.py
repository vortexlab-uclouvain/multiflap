import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import aero_package.bird_model as bm

class RetrievingAero:
    """Retrieve the aerodynamic variables of the limit cycle

    It is a fully post-processing routine which comes with multiflap and takes
    the periodic orbit as input, and re-run the aerodynamic modelaccordingly.

    Input:
        bird_obj: the bird object built with the constructor
        time_steps: the time array of the limit cycle.



    """

    def __init__(self, bird_obj, time_steps):
        self.bird_obj = bird_obj
        self.time_steps = time_steps

    def postrun_aeroforces(self, periodic_orbit):
        """Re-run the aerodynamic model

        """
        timeArray = self.time_steps
        time_steps = len(timeArray)
        # =============================================================================
        #     Initialization of forces array
        # =============================================================================
        Fx = np.zeros(time_steps)
        Fy = np.zeros(time_steps)
        Fz = np.zeros(time_steps)
        My = np.zeros(time_steps)
        F_tail = np.zeros(time_steps)
        M_wing = np.zeros(time_steps)
        M_tail = np.zeros(time_steps)
        M_drag = np.zeros(time_steps)
        M_lift = np.zeros(time_steps)
        P_ind = np.zeros(time_steps)
        for i in range (time_steps):
            [_, Fy[i], Fz[i], My[i], F_tail[i], _, _, _, _, P_ind[i]] = self.bird_obj.get_aeroforces(periodic_orbit[i], timeArray[i])
        return Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift, P_ind

    def parasitic_power(self, periodic_orbit):
        timeArray = self.time_steps
        time_steps = len(timeArray)
        P_par = np.zeros(len(self.time_steps))
        for i in range (time_steps):
            u = periodic_orbit[i, 0]
            w = periodic_orbit[i, 1]
            norm_v = np.sqrt(u**2 + w**2)
            cd_b = self.get_body_coeff(norm_v)
            s_b = 0.00813*self.bird_obj.mass**(2/3)
            drag_body = (0.5*self.bird_obj.rho*(cd_b*s_b)*(norm_v)**2)
            P_par[i] = drag_body*norm_v
        return P_par

    def get_body_coeff(self, v):
        rho = 1.225
        mu = 1.81*10**(-5)
          # Body surface
        S_b = 0.00813*self.bird_obj.mass**(2/3)

        # Diameter of the body
        d = np.sqrt(4*S_b/(np.pi))

        FR_t = 6.0799*self.bird_obj.mass**(0.1523) #fineness ratio l/d (only trunk of the bird)

        # Reynolds nuself.massber
        Re = (rho*d*v)/mu

        # MayBury Thesis
        CD_b = 66.6*self.bird_obj.mass**(-0.511)*FR_t**(0.915)*S_b**(1.063)*Re**(-0.197)

        return CD_b
