'''
file postrunning.py
@author Gianmarco Ducci
@copyright Copyright © UCLouvain 2020

multiflap is a Python tool for finding periodic orbits and assess their stability via the Floquet multipliers.

Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>

List of the contributors to the development of multiflap, Description and complete License: see LICENSE and NOTICE files.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
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
        Drag_tail = np.zeros(time_steps)
        M_wing = np.zeros(time_steps)
        M_tail = np.zeros(time_steps)
        M_drag = np.zeros(time_steps)
        M_lift = np.zeros(time_steps)
        P_ind = np.zeros(time_steps)
        F_pro = np.zeros(time_steps)
        P_pro = np.zeros(time_steps)
        angle_of_attack = []
        print(self.bird_obj.tail_opening)
        for i in range (time_steps):
            [_, Fy[i], Fz[i], My[i], F_tail[i],Drag_tail[i], M_wing[i], M_tail[i], _, _, P_ind[i], F_pro[i], aoa] = self.bird_obj.get_aeroforces(periodic_orbit[i], timeArray[i])
            aoa = np.rad2deg(aoa)
            P_pro[i] = F_pro[i]*periodic_orbit[i][0]
            angle_of_attack.append(aoa)
        return Fx, Fy, Fz, My, F_tail, Drag_tail, M_wing, M_tail, M_drag, M_lift, P_ind, P_pro, angle_of_attack
