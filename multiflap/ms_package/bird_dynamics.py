import numpy as np  # Import NumPy
from aero_package.flapping_forces import get_flappingforces
import scipy.integrate as ode
import time
from math import sin, cos
dim = 4

def dynamics(self, x0, t):
        """
        ODE system for the Equation of Motion on the longitudinal plane.
        This function will be passed to the numerical integrator

        Inputs:
            x0: initial values [u, w, q, theta]
            t: time

        Outputs:
            x_dot: velocity vector
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
            x0: initial values [u, w, q, theta]
            t: time

        Outputs:
            aereodynamic_forces : components of the aereod. forces and moments
        """
        dt = 0.25*1e-4
        wing_state = self.get_wingstate(t)
        wing_state_dt = self.get_wingstate(t-dt)

        [leadingedge,
         trailingedge] = self.get_wingenvelope(wing_state)

        [leadingedge_dt,
         trailingedge_dt] = self.get_wingenvelope(wing_state_dt)

        lifting_line =  self.get_liftingline(leadingedge,
                                             trailingedge)

        lifting_line_dt =  self.get_liftingline(leadingedge_dt,
                                                trailingedge_dt)

        [line, chordir, updir, chord] = self.merge_lines(lifting_line)

        [line_dt, _, _, _] = self.merge_lines(lifting_line_dt)

        v_kinematics = get_kinematicvelocity(line, line_dt, dt)

        # call the aereodynamic solver for the lifting line get_flappingforces
        [Fx,
         Fy,
         Fz,
         My,
         F_tail,
         M_wing,
         M_tail,
         M_drag,
         M_lift] = get_flappingforces(x0, v_kinematics,
                                      line, chordir, updir, chord)

        aereodynamic_forces = [Fx, Fy, Fz, My,
                               F_tail, M_wing, M_tail, M_drag, M_lift]
        return aereodynamic_forces


def get_kinematicvelocity(line, line_dt, dt):

       line_c = line[:, 1:-1:2]
       line_c2 = line_dt[:, 1:-1:2]
       v_kin =  (line_c - line_c2)/dt

       return v_kin

def get_stability_matrix(self, x0, t):

    """
    Stability matrix of the ODE system for the longitudinal plane

    Inputs:
        x0: initial condition [u_0, w_0, q_0, theta_0]
    Outputs:
        A: Stability matrix evaluated at x0. (dxd) dimension
        A[i, j] = dv[i]/dx[j]
    """
    m = self.mass
    g = self.g
    perturbation = 1e-3
    x0_u = x0 + [x0[0]*perturbation, 0, 0, 0]
    x0_w = x0 + [0, x0[1]*perturbation, 0, 0]
    x0_q = x0 + [0, 0, x0[2]*perturbation, 0]
    [_, Fy, Fz, My, F_tail, _, _, _, _] = self.get_aeroforces(x0, t)
    # force evaluation, perturbation along 'u'
    [_, Fyu, Fzu, Myu, F_tailu, _, _, _, _] = self.get_aeroforces(x0_u, t)
    # force evaluation, perturbation along 'w'
    [_, Fyw, Fzw, Myw, F_tailw, _, _, _, _] = self.get_aeroforces(x0_w, t)
    # force evaluation, perturbation along 'q'
    [_, Fyq, Fzq, Myq, F_tailq, _, _, _, _] = self.get_aeroforces(x0_q, t)


    # Derivatives of Fz with respect to the state space variables
    dFzu_dU = (Fzu - Fz)/(x0[0]*perturbation)
    dFzw_dW = (Fzw - Fz)/(x0[1]*perturbation)
    dFzq_dq = (Fzq - Fz)/(x0[2]*perturbation)

    # Derivatives of Fy with respect to the state space variables
    dFyu_dU = (Fyu - Fy)/(x0[0]*perturbation)
    dFyw_dW = (Fyw - Fy)/(x0[1]*perturbation)
    dFyq_dq = (Fyq - Fy)/(x0[2]*perturbation)

    # Derivatives of F_tail with respect to the state space variables
    dFytail_du = (F_tailu - F_tail)/(x0[0]*perturbation)
    dFytail_dw = (F_tailw - F_tail)/(x0[1]*perturbation)
    dFytail_dq = (F_tailq - F_tail)/(x0[2]*perturbation)

    # Derivatives of My with respect to the state space variables
    dMy_du = (Myu - My)/(x0[0]*perturbation)
    dMy_dw = (Myw - My)/(x0[1]*perturbation)
    dMy_dq = (Myq - My)/(x0[2]*perturbation)
    u, w, q, theta = x0

    A = np.array([[-dFzu_dU/m, -q - dFzw_dW/m,
                   -w - dFzq_dq/m, -g*np.cos(theta)],
                  [q - dFyu_dU/m - dFytail_du/m, -dFyw_dW/m - dFytail_dw/m,
                   u - dFyq_dq/m - dFytail_dq/m, -g*np.sin(theta)],
                  [dMy_du/0.1, dMy_dw/0.1, dMy_dq/0.1, 0],
                  [0, 0, 1, 0]], float)
    return A
