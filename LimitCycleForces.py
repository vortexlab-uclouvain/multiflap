import numpy as np
from FlappingForcesDev import FlappingForces
import settings as settings
import matplotlib.pyplot as plt
import matplotlib as mpl

def ForceRetrieving_CaseName(case_name, **kinematics):
    results_directory = '/Users/gducci/UCL/PROJECT/Simulations/NoIteration/'+ case_name+'/Results'
    periodic_orbit_filename = results_directory+'/complete_solution.npy'
    periodic_orbit = np.load(periodic_orbit_filename)
    periodic_orbit = periodic_orbit.reshape(-1, periodic_orbit.shape[2])
    time_steps = np.size(periodic_orbit, 0)
    timeArray = np.linspace(0, (1/settings.frequency), time_steps)
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
    u = periodic_orbit[:,0]
    w = periodic_orbit[:,1]
    q = periodic_orbit[:,2]
    theta = periodic_orbit[:,3]
    for i in range (time_steps):
        Fx[i], Fy[i], Fz[i], My[i], F_tail[i], M_wing[i], M_tail[i], M_drag[i], M_lift[i] = FlappingForces(timeArray[i],u[i],w[i],q[i],theta[i], **kinematics)
    return Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift, timeArray

def ForceRetrieving(periodic_orbit, **kinematics):
    periodic_orbit = periodic_orbit.reshape(-1, periodic_orbit.shape[2])
    time_steps = np.size(periodic_orbit, 0)
    timeArray = np.linspace(0, (1/settings.frequency), time_steps)
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
    u = periodic_orbit[:,0]
    w = periodic_orbit[:,1]
    q = periodic_orbit[:,2]
    theta = periodic_orbit[:,3]
    for i in range (time_steps):
        Fx[i], Fy[i], Fz[i], My[i], F_tail[i], M_wing[i], M_tail[i], M_drag[i], M_lift[i] = FlappingForces(timeArray[i],u[i],w[i],q[i],theta[i], **kinematics)
    return Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift



if __name__ == "__main__":
    case_name = 'Sim_00_v1'
    
    [Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift, timeArray] = ForceRetrieving_CaseName(case_name, amp_shoulder_z = np.deg2rad(42), amp_shoulder_y = np.deg2rad(20), off_shoulder_y = -np.deg2rad(15), tail_opening=0)
    fig1 = plt.figure()
    ax1 = fig1.gca() 
    fig1.suptitle('Moments', fontsize=18)
    plt.xlabel('1/T', fontsize=14)
    plt.ylabel('$M_i(t)$', fontsize=14)
#    ax1.plot(timeArray/0.25, Fy/(9.81*1.2), '-', label = "Moment wing")
##    ax1.plot(timeArray, Moment_tail, '-', label = "Moment tail")
    ax1.plot(timeArray, Moment_total, '-', label = "Moment")
    
    ax1.legend()
#    ax1.set_xlim(0, timeArray[-1]/0.25)
    ax1.grid(True)
    
    
    fig1 = plt.figure()
    ax1 = fig1.gca() 
    fig1.suptitle('Moments', fontsize=18)
    plt.xlabel('1/T', fontsize=14)
    plt.ylabel('$M_i(t)$', fontsize=14)
    ax1.plot(timeArray, Moment_drag, '-', label = "Moment drag")
    ax1.plot(timeArray, Moment_lift, '-', label = "Moment lift")
    ax1.plot(timeArray, Moment_tail, '-', label = "Moment tail")
    
    ax1.legend()
    ax1.set_xlim(0, timeArray[-1])
    ax1.grid(True)

    