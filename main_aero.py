import numpy as np  # Import NumPy
from flapping_forces import get_flappingforces
import settings as settings
import scipy.integrate as ode
import time
from math import sin, cos
from settings import SimulationsSettings
import bird_model as bm
from bird import Shoulder, Elbow, Wrist, Joint


def get_aeroforces(self, x0, t):
        """
        Returns the aereodynamic forces at a time t
        1. Get the wing state (angles of each joint at time t)
        2. Based on the joint angles, exctracts the wing envelope
        3. From the envelope, exctracts the lifing line
        4. The lifting line is mirrored for the symmetric wing
        5. Lifting line properties and bird velocities are
            passed to the liftin line solver

        Inputs:
        x0: initial values (that corresponds to the flow velocity)
        x0 = (u, w, q, theta)
        t: time

        Outputs:
        aereodynamic_forces : vector of the components iof the aereodynamic_forces and moments
        """
        dt = 0.25*1e-4
        wing_state = self.get_wingstate(t)
        wing_state_dt = self.get_wingstate(t-dt)
        [leadingedge, trailingedge] = self.get_wingenvelope(wing_state)
        [leadingedge_dt, trailingedge_dt] = self.get_wingenvelope(wing_state_dt)
        lifting_line =  self.get_liftingline(leadingedge, trailingedge)
        lifting_line_dt =  self.get_liftingline(leadingedge_dt, trailingedge_dt)
        [line, chordir, updir, chord] = self.merge_lines(lifting_line)
        [line_dt, _, _, _] = self.merge_lines(lifting_line_dt)
        v_kinematics = get_kinematicvelocity(line, line_dt, dt)

        # call the aereodynamic solver for the lifting line get_flappingforces
        [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = get_flappingforces(x0, v_kinematics, line, chordir, updir, chord)
        aereodynamic_forces = Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift
        return aereodynamic_forces


def get_kinematicvelocity(line, line_dt, dt):

       line_c = line[:, 1:-1:2]
       line_c2 = line_dt[:, 1:-1:2]
       v_kin =  (line_c - line_c2)/dt

       return v_kin

# generate bird kinematics by calling "bird" module
bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                         axis_y=Joint(-np.deg2rad(19),np.deg2rad(20),np.pi/2),
                         axis_z=Joint(0,np.deg2rad(42),np.pi))
bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                   axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                   axis_z=Joint(0.,0.,0.))

mybird = bm.BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist)


time_initial = 0.
time_final = 0.25

time_array = np.linspace(time_initial, time_final, 100)
x0 = [16, 0, 0, 0]
lift = []
for i in range(len(time_array)):

    mysolution = mybird.get_aeroforces(x0, time_array[i])
    lift.append(mysolution[1])


lift = np.asanyarray(lift)

import matplotlib.pyplot as plt

plt.plot(time_array, lift)
plt.show()
