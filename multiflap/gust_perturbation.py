import numpy as np  # Import NumPy
import RungeKutta as rk
from FlappingForcesDev import FlappingForces
import settings as settings
import scipy.integrate as ode
from LimitCycleForces import ForceRetrieving
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from scipy import interpolate

dim = 4

f = settings.frequency

count = 0
# Parameters:
g = 9.81
mass = 1.2

# =============================================================================
# Initialization of the aerodynamic forces to be extracted
# =============================================================================

gust_lift = []
gust_drag = []
gust_moment = []
gust_moment_lift = []
gust_moment_drag = []
gust_moment_tail = []
gust_moment_wing = []
# =============================================================================
# Equations of motion
# =============================================================================


def Velocity(ssp, t, **kinematics):     # Standard EoM
    # Read inputs:
    u, w, q, theta = ssp  # Read state space points
    u = ssp[0]
    w = ssp[1]
    q = ssp[2]
    theta = ssp[3]
    # Longitudinal equation of motion:
    [Fx, Fy, Fz, My, F_tail, M_wing,
     M_tail, M_drag, M_lift] = FlappingForces(t, u, w, q, theta, **kinematics)
    dudt = -q*w - g*np.sin(theta) - Fz/mass
    dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
    dqdt = My/0.1
    dthetadt = q
    # Collect Equations of motion in a single NumPy array:

    vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector

    return vel


def birdEqn_py_Gust(t, ssp, **kinematics):
        """
        State space velocity function for the Equation of Motion of the longitudinal plane

        Inputs:
        ssp: State space vector
        ssp = (u, w, q, theta)
        t: Time

        Outputs:
        vel: Time derivative of ssp.
        """
        #Parameters:
        g = 9.81
        mass = 1.2
        #Read inputs:
        u, w, q, theta  = ssp  # Read state space points
        u = ssp[0]
        w = ssp[1]
        q = ssp[2]
        theta = ssp[3]
        w_gust = gaussian_gust(t)

        # Longitudinal equation of motion:
        [Fx, Fy, Fz, My, F_tail, M_wing, M_tail,
         M_drag, M_lift] = FlappingForces(t, u, w+w_gust, q,
                                         theta, **kinematics)
        dudt = -q*w - g*np.sin(theta) - Fz/mass
        dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
        dqdt =  My/0.1
        dthetadt = q                # Collect Equations of motion in a single NumPy array:

        vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector            
        return vel

def Velocity_Gust(ssp, t, **kinematics):
    # EoM with perturbation (gust)
    u, w, q, theta = ssp  # Read state space points
    u = ssp[0]
    w = ssp[1]
    q = ssp[2]
    theta = ssp[3]
    w_gust = gaussian_gust(t)
    [Fx, Fy, Fz, My, F_tail, M_wing, M_tail,
     M_drag, M_lift] = FlappingForces(t, u, w+w_gust, q,
                                      theta, **kinematics)
    gust_lift.append(Fy)
    gust_drag.append(Fz)
    gust_moment.append(My)
    gust_moment_lift.append(M_lift)
    gust_moment_drag.append(M_drag)
    gust_moment_tail.append(M_tail)
    gust_moment_wing.append(M_wing)
    # Longitudinal equation of motion:

    dudt = -q*w - g*np.sin(theta) - Fz/mass
    dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
    dqdt = My/0.1
    dthetadt = q
    # Collect Equations of motion in a single NumPy array:

    vel = np.array([dudt, dwdt, dqdt, dthetadt], float)
    return vel


def gaussian_gust(t):
    sigma = .05
    t_0 = .3  # sec, max of the gaussian
    a_0 = 2  # peak of the gaussian (max(w(t)))
    w = a_0*np.e**(-0.5*(((t-t_0)/sigma)**2))
    return w


if __name__ == "__main__":
    settings.amplitude_shoulder_z = np.deg2rad(39.45047827064355)
    settings.offset_shoulder_y = -np.deg2rad(26)
    settings.tail_opening = np.deg2rad(40)
    force_retrieving = True
#    case_name = 'TestCase12b'

    if force_retrieving is True:
        case_name = 'NonLevel'
        results_directory = '/Users/gducci/UCL/PROJECT/Simulations/ResultsTail/LevelSimulations/TailOpening_40/SweepAmplitude_20/SweepOff_Neg-26'

        periodic_orbit_filename = results_directory+'/complete_solution.npy'
        periodic_orbit = np.load(periodic_orbit_filename)
        u_0 = periodic_orbit[0, 0][0]  # U-velocity initial condition
        w_0 = periodic_orbit[0, 0][1]  # W-velocity initial condition
        q_0 = periodic_orbit[0, 0][2]  # Q-velocity initial condition
        theta_0 = periodic_orbit[0, 0][3]  # Theta-angle initial condition
        ssp0 = np.array([u_0, w_0, q_0, theta_0], float)
        periodic_orbit = periodic_orbit.reshape(-1, periodic_orbit.shape[2])
        complete_solution = np.load(periodic_orbit_filename)
# ====================================================================
# Load forces of the limit cycle solution
# ====================================================================
        lift_filename = results_directory+'/Lift_coupled_v2.npy'
        force_tail_filename = results_directory+'/Force_tail.npy'
        drag_filename = results_directory+'/Drag_coupled_v2.npy'
        Moment_total_filename = results_directory+'/Moment_total.npy'
        Moment_wing_filename = results_directory+'/Moment_wing.npy'
        Moment_tail_filename = results_directory+'/Moment_tail.npy'
        Moment_drag_filename = results_directory+'/Moment_drag.npy'
        Moment_lift_filename = results_directory+'/Moment_lift.npy'
        lift_coupled = np.load(lift_filename)
        drag_coupled = np.load(drag_filename)
        force_tail = np.load(force_tail_filename)
        Moment_total = np.load(Moment_total_filename)
        Moment_wing = np.load(Moment_wing_filename)
        Moment_tail = np.load(Moment_tail_filename)
        Moment_drag = np.load(Moment_drag_filename)
        Moment_lift = np.load(Moment_lift_filename)

        u_ref = periodic_orbit[:, 0]
        w_ref = periodic_orbit[:, 1]
        q_ref = periodic_orbit[:, 2]
        theta_ref = periodic_orbit[:, 3]
        eigenvalues_dir = results_directory+'/outfile_JacobianEigenvalues.npy'
        eigenvectors_dir = results_directory+'/outfile_JacobianEigenvector.npy'
        eigenvalues = np.load(eigenvalues_dir)
        eigenvectors = np.load(eigenvectors_dir)
    else:
        u_0 = 16.   # U-velocity initial condition
        w_0 = .0          # W-velocity initial condition
        q_0 = 0.  # Q-velocity initial condition
        theta_0 = 0.  # Theta-angle initial condition
        ssp0 = np.array([u_0, w_0, q_0, theta_0], float)

# ====================================================================
# Definition of the time array
# ====================================================================
    period_number = 5    # Number of periods over which integrate the EoM
    tInitial = 0                      # Initial time
    tFinal = period_number*(1/f)      # Final time (This is one period)
    Nt = 40*period_number                   # Discretisation of time array
    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array

# ====================================================================
# Repeating solution array, for the numbers of periods selected
# ====================================================================
    u_ref_solution = np.tile(u_ref[0:-1], period_number)
    u_ref_solution = np.append(u_ref_solution, u_ref[0])

    w_ref_solution = np.tile(w_ref[0:-1], period_number)
    w_ref_solution = np.append(w_ref_solution, w_ref[0])

    q_ref_solution = np.tile(q_ref[0:-1], period_number)
    q_ref_solution = np.append(q_ref_solution, q_ref[0])

    theta_ref_solution = np.tile(theta_ref[0:-1], period_number)
    theta_ref_solution = np.append(theta_ref_solution, theta_ref[0])

    theta_ref_solution = np.tile(theta_ref[0:-1], period_number)
    theta_ref_solution = np.append(theta_ref_solution, theta_ref[0])

    lift_coupled_ref = np.tile(lift_coupled[0:-1], period_number)
    lift_coupled_ref = np.append(lift_coupled_ref, lift_coupled[0])

    lift_tail_ref = np.tile(force_tail[0:-1], period_number)
    lift_tail_ref = np.append(lift_tail_ref, force_tail[0])

    moment_total_ref = np.tile(Moment_total[0:-1], period_number)
    moment_total_ref = np.append(moment_total_ref, Moment_total[0])

    moment_wing_ref = np.tile(Moment_wing[0:-1], period_number)
    moment_wing_ref = np.append(moment_wing_ref, Moment_wing[0])

    moment_tail_ref = np.tile(Moment_tail[0:-1], period_number)
    moment_tail_ref = np.append(moment_tail_ref, Moment_tail[0])

    tArray_ref = np.linspace(0, (1/f)*period_number, len(u_ref_solution))

    sspSolution = rk.RK2(Velocity, ssp0,tArray)
    sspSolution_gust = rk.RK2(Velocity_Gust, ssp0, tArray)
#    sspSolution_ivp = ode.solve_ivp(birdEqn_py_Gust, [tInitial, tFinal], ssp0, max_step=0.001)

# ====================================================================
# Gust solution
# ====================================================================
    velocity_U = sspSolution_gust[:, 0]
    velocity_W = sspSolution_gust[:, 1]
    Theta = sspSolution_gust[:, 3]

    lift = np.asanyarray(gust_lift)
    gust_moment = np.asanyarray(gust_moment)
    gust_moment_tail = np.asanyarray(gust_moment_tail)
    gust_moment_wing = np.asanyarray(gust_moment_wing)

    lift_2 = np.zeros(len(tArray))
    gust_moment_2 = np.zeros(len(tArray))
    gust_moment_tail_2 = np.zeros(len(tArray))
    gust_moment_wing_2 = np.zeros(len(tArray))
    lift_2[0:-1] = lift[::2]
    lift_2[-1] = lift[-1]

    gust_moment_2[0:-1] = gust_moment[::2]
    gust_moment_2[-1] = gust_moment[-1]

    gust_moment_tail_2[0:-1] = gust_moment_tail[::2]
    gust_moment_tail_2[-1] = gust_moment_tail[-1]

    gust_moment_wing_2[0:-1] = gust_moment_wing[::2]
    gust_moment_wing_2[-1] = gust_moment_wing[-1]
    u_hor = velocity_U*np.cos(Theta) + velocity_W*np.sin(Theta)
    v_vert = velocity_U*np.sin(Theta) - velocity_W*np.cos(Theta)
    z_initial = 10.
    x = ((np.cumsum(u_hor))/len(u_hor)/4)
    z = z_initial + ((np.cumsum(v_vert))/len(u_hor)/4)

    u_hor_ref = sspSolution[:,0]*np.cos(sspSolution[:,3]) + sspSolution[:,1]*np.sin(sspSolution[:,3])
    v_vert_ref = sspSolution[:,0]*np.sin(sspSolution[:,3]) - sspSolution[:,1]*np.cos(sspSolution[:,3])
    z_initial_ref = 10.
    x_ref = ((np.cumsum(u_hor_ref))/len(u_hor_ref)/4)
    z_ref = z_initial_ref + ((np.cumsum(v_vert_ref))/len(u_hor_ref)/4)
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import proj3d

#    import matplotlib.animation as animation

    fig1 = plt.figure()
    plt.plot(x,  z)
    plt.plot(x_ref,  z_ref)
#    plt.plot(tArray,  gust_moment_wing_2)

#    plt.plot(tArray,  velocity_U, color='red')
    plt.show()

    time_plot = np.linspace(0, 1, 1000)
    gust_plot = [gaussian_gust(time_plot[i]) for i in range(len(time_plot))]
    gust_plot = np.asarray(gust_plot)

    fig = plt.figure()
    ax_ref = fig.gca()
    ax_ref.plot(time_plot/0.25, gust_plot, linewidth=2.5, color='red')
    ax_ref.set_xlim(0, 1/0.25)
#    plt.plot(tArray,  lift_2, color='red')
    ax_ref.set_xlabel("$t/T$", fontsize=16)
    ax_ref.set_ylabel("$w_{gust}[m/s]$", fontsize=16)
    ax_ref.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.5))
    ax_ref.yaxis.set_major_locator(mpl.ticker.MultipleLocator(.5))
    plt.tick_params(labelsize=12)
    plt.grid()
    plt.show()
    # plt.savefig('/Users/gducci/UCL/MyWrittenDocs/
    # Report/ResultsTail/figures/'+
    # 'gust_gaussian.eps', format = 'eps')

    u_ref_mean = [np.mean(u_ref_solution)]*len(velocity_U)
    w_ref_mean = [np.mean(w_ref_solution)]*len(velocity_U)
    q_ref_mean = [np.mean(q_ref_solution)]*len(velocity_U)
    theta_ref_mean = [np.mean(theta_ref_solution)]*len(velocity_U)

    fig = plt.figure(figsize=(18, 16))
    gs = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=[1., 1., 1.,  1.])
    ax1 = fig.add_subplot(gs[0, 0])
#    ax1.grid()
    ax1.set_xlim(0, tArray[-1])
    ax1.plot(tArray, velocity_U, color='blue')
    ax1.plot(tArray_ref, u_ref_solution, color='red', alpha=0.3)
    ax1.plot(tArray, u_ref_mean, color='red')
    ax1.set_ylabel('$u[m/s]$', fontsize=16)
    ax1.grid()
    plt.tick_params(labelsize=12)

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_xlim(0, tArray[-1])
    ax2.set_ylabel('$w[m/s]$', fontsize=16)
    ax2.plot(tArray, velocity_W, color='blue')
    ax2.plot(tArray_ref, w_ref_solution, color='red', alpha=0.3)
    ax2.plot(tArray, w_ref_mean, color='red')
    ax2.grid()
    plt.tick_params(labelsize=12)

    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_xlim(0, tArray[-1])
    ax3.grid()
    ax3.plot(tArray, sspSolution_gust[:, 2], color='blue')
    ax3.plot(tArray_ref, q_ref_solution, color='red', alpha=0.3)
    ax3.plot(tArray, q_ref_mean, color='red')

    ax3.set_ylabel('$q[rad/s]$', fontsize=16)
    plt.tick_params(labelsize=12)

    ax_theta = fig.add_subplot(gs[3, 0])
    ax_theta.set_xlim(0, tArray[-1])
    ax_theta.grid()
    ax_theta.plot(tArray, sspSolution_gust[:, 3], color='blue')
    ax_theta.plot(tArray_ref, theta_ref_solution, color='red', alpha=0.3)
    ax_theta.plot(tArray, theta_ref_mean, color='red')

    ax_theta.set_ylabel('$\Theta[rad]$', fontsize=16)
    ax_theta.set_xlabel("$t[s]$", fontsize=16)
    plt.tick_params(labelsize=12)
#    plt.savefig('/Users/gducci/UCL/MyWrittenDocs/Report/ResultsTail/figures/'+'gust_longterm.pdf', format = 'pdf')

    fig = plt.figure(figsize=(18, 16))
    gs = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[1., 1., 1.])
    ax_wing = fig.add_subplot(gs[0, 0])
#    ax1.grid()
    ax_wing.set_xlim(0, tArray[-1])
    ax_wing.plot(tArray, gust_moment_wing_2, color='blue')
    ax_wing.plot(tArray_ref, moment_wing_ref, color='red', alpha=0.3)
    ax_wing.set_ylabel('$M_{wing}[Nm]$', fontsize=16)
    ax_wing.grid()
    plt.tick_params(labelsize=12)

    ax_tail = fig.add_subplot(gs[1, 0])
    ax_tail.set_xlim(0, tArray[-1])
    ax_tail.set_ylabel('$M_{tail}[Nm]$', fontsize=16)
    ax_tail.plot(tArray, gust_moment_tail_2, color='blue')
    ax_tail.plot(tArray_ref, moment_tail_ref, color='red', alpha=0.3)
    ax_tail.grid()
    ax_tail.tick_params(labelsize=12)

    ax_total = fig.add_subplot(gs[2, 0])
    ax_total.set_xlim(0, tArray[-1])
    ax_total.grid()
    ax_total.plot(tArray, gust_moment_2, color='blue')
    ax_total.plot(tArray_ref, moment_total_ref, color='red', alpha=0.3)

    ax_total.set_ylabel('$M_{total}[Nm]$', fontsize=16)
    plt.tick_params(labelsize=12)
    ax_total.set_xlabel("$t[s]$", fontsize=16)
    plt.tick_params(labelsize=12)
    plt.show()
    # plt.savefig('/Users/gducci/UCL/MyWrittenDocs/Report/ResultsTail/figures/'+'gust_longterm_moments.pdf', format = 'pdf')

    lines = []

    gust_vel = [gaussian_gust(tArray[i]) for i in range(len(tArray))]
    gust_vel = np.asarray(gust_vel)

# ===========================================================
# Making the two solutions of equal dimension 
# ===========================================================
    def interpolator(t, array):
        f = interpolate.interp1d(tArray_ref[:], array[:], kind='cubic')
        var_ref = f(t)
        return var_ref

    u_ref_shaped = [interpolator(tArray[i], u_ref_solution)
                    for i in range(len(tArray))]
    w_ref_shaped = [interpolator(tArray[i], w_ref_solution)
                    for i in range(len(tArray))]
    q_ref_shaped = [interpolator(tArray[i], q_ref_solution)
                    for i in range(len(tArray))]
    theta_ref_shaped = [interpolator(tArray[i], theta_ref_solution)
                        for i in range(len(tArray))]
    moment_tail_reshaped = [interpolator(tArray[i], moment_tail_ref)
                            for i in range(len(tArray))]

#    for i in range(len(tArray)):
#        print(i)
#        head = i
#        line1,  = ax1.plot(tArray[:i], w_ref_shaped[:i],
#                           color='blue')
#        line1_inst, = ax1.plot(tArray[:i], velocity_W[:i],
#                             color='red')
#
#        line1_head, = ax1.plot(tArray[i], w_ref_shaped[i],
#                             color='blue', marker='o', markeredgecolor='blue')
#        line1_inst_head, = ax1.plot(tArray[i], velocity_W[i],
#                         color='red', marker='o', markeredgecolor='red')
#
#        line2,  = ax2.plot(tArray[:i], u_ref_shaped[:i],
#                        color='blue')
#        line2_inst, = ax2.plot(tArray[:i], velocity_U[:i],
#                        color='red')
#        line2_head_inst_0, = ax2.plot(tArray[i], u_ref_shaped[i],
#                            color='blue', marker='o', markeredgecolor='blue')
#
#        line2_head, = ax2.plot(tArray[i], velocity_U[i],
#                           color='blue', marker='o', markeredgecolor='blue')
#        line2_head_inst, = ax2.plot(tArray[i], velocity_U[i],
#                           color='red', marker='o', markeredgecolor='red')
#
#
#        line3,  = ax3.plot(tArray[:i], gust_vel[:i], color='red')
#        line3_head_inst,  = ax3.plot(tArray[i], gust_vel[i],
#                           color='red', marker='o', markeredgecolor='red')
#
#        lines.append([line1, line1_inst, line1_head, line1_inst_head, line2,
#                      line2_inst,line2_head,line2_head_inst, line2_head_inst_0,
#                      line3, line3_head_inst])
#
#    plt.tight_layout()
#    ani = animation.ArtistAnimation(fig, lines, interval=50, blit=True)
#    Writer = animation.writers['ffmpeg']
#    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#    ani.save('/Users/gducci/Desktop/gust_tail.mp4')
#
#    fig = plt.figure()
#    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1.,1.])
#    ax4 = fig.add_subplot(gs[0, 0])
#    ax4.grid()
#    ax4.set_ylabel('$M_{wing}, M_{tail}[Nm]$')
#    ax5 = fig.add_subplot(gs[1, 0])
#    ax5.set_ylabel('$M_{total}[Nm]$')
#    ax5.grid()
#
#    lines = []
#
#    for i in range(len(tArray[0:300])):
#        print(i)
#        head = i
#        line4,  = ax4.plot(tArray[:i], gust_moment_tail_2[:i],
#                           color='blue')
#        line4_inst, = ax4.plot(tArray[:i], gust_moment_wing_2[:i],
#                           color='red')
#        line4_inst_tail, = ax4.plot(tArray[:i], moment_tail_reshaped[:i],
#                           color='blue', alpha=0.3)
#
#        line4_head, = ax4.plot(tArray[i], gust_moment_tail_2[i],
#                           color='blue', marker='o', markeredgecolor='blue')
#        line4_inst_head,= ax4.plot(tArray[i], gust_moment_wing_2[i],
#                       color='red', marker='o', markeredgecolor='red')
#
#
#        line5,  = ax5.plot(tArray[:i], gust_moment_2[:i],
#                       color='blue')
#        line5_head_inst_0, = ax5.plot(tArray[i], gust_moment_2[i],
#                           color='blue', marker='o', markeredgecolor='blue')
#
#
#
#
#
#        lines.append([line4,line4_inst, line4_head,line4_inst_tail, line4_inst_head, line5, line5_head_inst_0])
#
#
#
#
#    plt.tight_layout()
#    ani_moment = animation.ArtistAnimation(fig, lines, interval=50, blit=True)

# =============================================================================
#     Making reference solution and gust solution of equal dimension
# =============================================================================


#    plt.plot(x_ref,  z_ref)
