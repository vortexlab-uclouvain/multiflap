import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import aero_package.bird_model as bm

class RetrievingAero:

    def __init__(self, bird_obj, time_steps):
        self.bird_obj = bird_obj
        self.time_steps = time_steps

    def ForceRetrieving(self, periodic_orbit):
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
