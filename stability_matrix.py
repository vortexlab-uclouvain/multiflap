import numpy as np
from bird_dynamics import get_aeroforces

def get_stability_matrix(self, x0, t):

    """
    Stability matrix of the ODE system for the longitudinal plane

    Inputs:
    x0: initial condition [u_0, w_0, q_0, theta_0]
    Outputs:
    A: Stability matrix evaluated at x0. (dxd) dimension
       A[i, j] = dv[i]/dx[j]
    """
    g = 9.81
    m = 1.2
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

    A = np.array([[-dFzu_dU/m, -q - dFzw_dW/m, -w - dFzq_dq/m, -g*np.cos(theta)],
                  [q - dFyu_dU/m - dFytail_du/m, -dFyw_dW/m - dFytail_dw/m, u - dFyq_dq/m - dFytail_dq/m, -g*np.sin(theta)],
                  [dMy_du/0.1, dMy_dw/0.1, dMy_dq/0.1, 0],
                  [0, 0, 1, 0]], float)
    return A
