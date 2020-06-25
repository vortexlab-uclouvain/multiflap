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
        State space velocity function for the Equation of Motion of the longitudinal plane

        Inputs:
        ssp: State space vector
        ssp = (u, w, q, theta)
        t: Time

        Outputs:
        vel: Time derivative of ssp.
        """
        #Read inputs:
        u, w, q, theta  = x0  # Read state space points
        # Longitudinal equation of motion:

        af = self.get_aeroforces(x0, t)

        # bird body dynamics
        dudt = -self.g*w - self.g*np.sin(theta) - af.Fz/self.mass
        dwdt = q*u + self.g*np.cos(theta) - af.Fy/self.mass - af.F_tail/self.mass
        dqdt =  af.My/0.1
        dthetadt = q
        # Collect Equations of motion in a single NumPy array:

        vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector

        return vel

def get_aeroforces(self, x0, t):
        # aereodynamic forces
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
        [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = get_flappingforces(x0, v_kinematics, line, chordir, updir, chord)
        aereodynamic_forces = Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift
        return aereodynamic_forces


def get_kinematicvelocity(line, line_dt, dt):

       line_c = line[:, 1:-1:2]
       line_c2 = line_dt[:, 1:-1:2]
       v_kin =  (line_c - line_c2)/dt

       return v_kin



if __name__ == "__main__":
#    settings.amplitude_shoulder_z = np.deg2rad(39.45047827064355)
#    settings.offset_shoulder_y = -np.deg2rad(26)
#    settings.tail_opening = np.deg2rad(40)

    force_retrieving =  True
#    case_name = 'TestCase12b'

    if force_retrieving ==  True:
        case_name = 'NonLevel'
        results_directory = '/Users/gducci/UCL/PROJECT/Simulations/ResultsPaper/M10_ref/Results'

        periodic_orbit_filename = results_directory+'/complete_solution.npy'
        periodic_orbit = np.load(periodic_orbit_filename)
        u_0 =  periodic_orbit[0,0][0]    # U-velocity initial condition
        w_0 =  periodic_orbit[0,0][1]         # W-velocity initial condition
        q_0 =  periodic_orbit[0,0][2]    # Q-velocity initial condition
        theta_0 = periodic_orbit[0,0][3]     # Theta-angle initial condition
        ssp0 = np.array([u_0, w_0, q_0, theta_0], float) # Vector of initial conditions
#        periodic_orbit = periodic_orbit.reshape(-1, periodic_orbit.shape[2])

    else:
        u_0 =  16.   # U-velocity initial condition
        w_0 = .0          # W-velocity initial condition
        q_0 = 0. # Q-velocity initial condition
        theta_0 = 0.  # Theta-angle initial condition
        ssp0 = np.array([u_0, w_0, q_0, theta_0], float) # Vector of initial conditions


    # -------------------------------------------------------------------------
    #  Definition of the time array
    # -------------------------------------------------------------------------
    period_number = 1       # Number of periods over which integrate the EoM
    tInitial = 0                      # Initial time
    tFinal = period_number*(1/f)      # Final time (This is one period)
    Nt = 10*period_number                   # Discretisation of time array
    tArray = np.linspace(tInitial, tFinal, 80)  # Time array
    tail_op = np.deg2rad(0)

    sspSolution_V0 = rk.RK2(Velocity, ssp0,tArray)
    test = Flow(ssp0, 0, tFinal, 80)

#    Jac = JacobianNumerical(ssp0, tArray[40], 0.25, endpoint=True)
#
#    eignevalues_numerical, eigenvector_numerical = np.linalg.eig(Jac)
#    np.save(results_directory+'/eigenvalues_jacobian_num', eignevalues_numerical)
#    end_Jac = time.time()
#    time_Jac = end_Jac - start_Jac
#    print(Jacobian)
#    print("RK2 = ", time_Jac)
#



