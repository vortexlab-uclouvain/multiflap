import numpy as np  # Import NumPy
import RungeKutta as rk
from flapping_forces import get_flappingforces
import settings as settings
import scipy.integrate as ode
import time
from math import sin, cos
from settings import SimulationsSettings
dim = 4
f = SimulationsSettings().frequency

def dynamics(self, x0, t):
        """
        ODE system for the Equation of Motion on the longitudinal plane.
        This function will be passed to the numerical integrator

        Inputs:
        x0: initial values
        x0 = (u, w, q, theta)
        t: time

        Outputs:
        x_dot : velocity vector
        """
        #Read inputs:
        u, w, q, theta  = x0  # Read state space points

        # get the aereodynamic forces at time t
        [_, Fy, Fz, My, F_tail, _, _, _, _] = self.get_aeroforces(x0, t)

        # bird body dynamics
        dudt = -q*w - self.g*np.sin(theta) - Fz/self.mass
        dwdt = q*u + self.g*np.cos(theta) - Fy/self.mass - F_tail/self.mass
        dqdt =  My/0.1
        dthetadt = q

        # collecting x_dot components in a single array:
        x_dot = np.array([dudt, dwdt, dqdt, dthetadt], float)

        return x_dot

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

