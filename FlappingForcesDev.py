import numpy as np  # Import NumPy
from BirdLine import BirdLine
from MergeLines import MergeLines
import Circulation 
import settings as settings
import TailModels as TailModel

def FlappingForces(t, u, w, theta, q, **kinematics):
    cl_alpha = 2*np.pi
    rho = 1.225 # Air density
    
    U = np.array([0, w, u])
    
    #angular_velocity = np.array([q, 0, 0])
    
    dt = 1e-6
    t2 = t - dt
    wingframe_position = settings.wingframe_position
    # -------------------------------------------------------------------------
    # Initialization of the vector force. It must have the same dimension of the time vector
    # -------------------------------------------------------------------------
    Fx = 0
    Fy = 0
    Fz = 0
    My = 0
    M_wing = 0 
    M_drag = 0
    M_lift = 0
    M_tail = 0
    F_tail_tot = 0
    # -------------------------------------------------------------------------
    # For each time step t, the bird line is calculated, and therefore gamma and forces
    # -------------------------------------------------------------------------
    x_ep = 0.05
    
    """
    Evaluation of the geometry at time t
    """
    
    # -------------------------------------------------------------------------
    # Call bird line function, to extract the Lifting line of right and left wing
    # -------------------------------------------------------------------------
    
    [line_right, up_direction_right, chord_direction_right, chord_right, 
     line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(t, **kinematics)
    
    nmid = np.ceil(2*x_ep/dx)
    # -------------------------------------------------------------------------
    # Geometry at time t
    # -------------------------------------------------------------------------
    [line, chordir, updir, chord] = MergeLines(line_right, up_direction_right, 
    chord_direction_right, chord_right, line_left,up_direction_left, chord_direction_left, chord_left,nmid,x_ep)
    
    
    """
    Evaluation of the geometry at time t2 = t - dt
    """
    # -------------------------------------------------------------------------
    # Call bird line function, to extract the Lifting line of right and left wing
    # -------------------------------------------------------------------------
    
    [line_right, up_direction_right, chord_direction_right, chord_right, 
     line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(t2, **kinematics)
    
    # -------------------------------------------------------------------------
    #  Geometry at time t2 = t-dt
    # -------------------------------------------------------------------------
    
    [line2, _, _, _] = MergeLines(line_right, up_direction_right,
    chord_direction_right, chord_right, line_left,up_direction_left, chord_direction_left, chord_left,nmid,x_ep)
    


    line_c = line[:,1:-1:2] # "Center" points: at the center of each interval --> every other point, starting from the second until the one before the last one
    line_c2 = line2[:,1:-1:2] # Center points at t-dt
    updir = updir[:,1:-1:2] # updir, chordir and chord at the center points
    chordir = chordir[:,1:-1:2]
    chord = chord[1:-1:2]
    line = line[:,::2] # "mid" points: points at the junction between two segments --> every other point from the first one to the last one

    # -------------------------------------------------------------------------
    #  u_inf = free stream velocity, equal over each profile (points)
    # -------------------------------------------------------------------------

    u_inf = np.array([U for j in range(np.size(line_c[0]))]).T
    #tangential_wing = (line_c - line_c2)/dt 
    velocity = u_inf - (line_c - line_c2)/dt 
    velocity_profile = np.zeros_like(velocity)
    line_direction = np.zeros_like(updir)        
    norm_velocity = np.zeros_like(chord)        
    direction_velocity = np.zeros_like(velocity)        
    angle_of_attack = np.zeros_like(chord)
    
    for j in range(np.size(line_c[0])):
        line_direction[:,j] = np.cross(updir[:,j], chordir[:,j])
        velocity_profile[:,j] = velocity[:,j] - line_direction[:,j]*np.sum(line_direction[:,j]*velocity[:,j])
        norm_velocity[j] = np.sqrt(np.sum(velocity_profile[:,j]*velocity_profile[:,j]))
        direction_velocity[:,j] = velocity_profile[:,j]/norm_velocity[j] 
        angle_of_attack[j] = np.arctan2(np.sum(direction_velocity[:,j]*updir[:,j]), np.sum(direction_velocity[:,j]*chordir[:,j]))
    tol = 0.01
    max_iterations = 50
    ds = np.zeros(len(chord))
    s = np.zeros((np.size(line[0])))
#    print(np.mean(np.rad2deg(angle_of_attack)))
    # ------------------------------------------------------------------------- 
    # Evaluation of ds length by distance point to point 
    # ds is the distance between two consecutives discretised points 
    # that belong to the lifting line
    # -------------------------------------------------------------------------

    for j in range(np.size(line[0])-1):
        ds[j] = np.sqrt(np.sum((line[:,j+1] - line[:,j])**2))
    
    for j in range(1,np.size(line[0])):
        s[j] = s[j-1] + ds[j-1]
        
    gamma = np.zeros((1,np.size(ds)))
    error = tol + 1
    error_p = np.inf
    count = 0
    
    while (error > tol) and (count < max_iterations):
        
        count = count + 1
        
        [gamma_new, v_downwash] = Circulation.Circulation(line[0,:], line[1,:], line_c[0,:], line_c[1,:],chord, angle_of_attack, cl_alpha, norm_velocity, gamma)
        
        error = np.mean(np.abs(gamma - gamma_new))/np.mean(np.abs(gamma_new))
        
        if error > error_p and count > 0:
            print("Error increasing after iteration ",count)
            break
        
        error_p = error
        
        gamma = gamma_new           
    
    gamma_filament = np.c_[gamma[:,0], gamma[:,1:] - gamma[0:,:-1], -gamma[0,-1]]
    gamma = np.reshape(gamma, np.size(gamma))
    gamma_filament = np.reshape(gamma_filament, np.size(gamma_filament))
    local_velocity = velocity + v_downwash
    lever_arm = np.zeros_like(line_c)
    BS_radius = np.zeros_like(line_c)   
# =============================================================================
#     Induced velocity on the tail
# =============================================================================
    [tail_span, AR_tail, NP_tail] = TailModel.tail_geometry(settings.tail_length, **kinematics)
    V_ind = 0
    dl_induced = line[:,1:] - line[:,:-1]
    for i in range (np.size(gamma)):
        lever_arm[:, i] = line_c[:, i] + wingframe_position         
        BS_radius[:,i] = (NP_tail - lever_arm[:, i])
        cross_prod = (np.cross(dl_induced[:,i],BS_radius[:,i]))
        V_ind = V_ind + ((gamma[i]/(4*np.pi))*cross_prod)/(np.linalg.norm(BS_radius[:,i])**3)
#        V_ind = V_ind + Circulation.InducedVel_Segment(NP_tail,line[:,i]+wingframe_position, line[:,i+1]+wingframe_position, gamma[i], 1e-13)

    lever_arm1 = np.zeros_like(line)
    BS_radius1 = np.zeros_like(line)
    for i in range (len(line)):
        lever_arm1[:, i] = line[:, i] + wingframe_position
        BS_radius1[:,i] = (NP_tail - lever_arm1[:, i])
        V_ind = V_ind + Circulation.InducedVel_Filament(BS_radius1[:,i],line[:,i], U, gamma_filament[i])
    U_tail = U + V_ind
    M = np.zeros_like(gamma)
    alpha_tail = np.arctan2(U_tail[1],U_tail[2])
    
    for j in range(np.size(gamma)):
        F = rho*gamma[j]*np.cross(local_velocity[:,j] , line_direction[:,j])
        lever_arm[:, j] = line_c[:, j] + wingframe_position
        F_moment = np.array([F[0]*ds[j], F[1]*ds[j], F[2]*ds[j]])
        Drag_moment = np.array([0., 0., F[2]*ds[j]])
        Lift_moment = np.array([0., F[1]*ds[j], 0.])

        Fx  =   Fx + F[0]*ds[j]         # Component along wing span
        Fy  =   Fy + F[1]*ds[j]         # Component of Lift (perpendicular to the free stream)
        Fz  =   Fz + F[2]*ds[j]         # Component of drag (along free stream)
        M_drag_j = np.cross(lever_arm[:,j], Drag_moment)
        M_lift_j = np.cross(lever_arm[:,j], Lift_moment)
        M_wing_j = np.cross(lever_arm[:,j], F_moment)
        
        M_wing = M_wing + M_wing_j[0]
        M_drag = M_drag + M_drag_j[0]
        M_lift = M_lift + M_lift_j[0]
        My = My + M[0]
    
    F_tail_tot = TailModel.delta_tail(U_tail, tail_span)
    
    M_tail = np.cross(NP_tail, F_tail_tot)
    My = M_wing + M_tail[0]
    return Fx, Fy, Fz, My, F_tail_tot[1], M_wing, M_tail[0], M_drag, M_lift



"""
The following allows to run the lifting line model without coupling it with the 
dynamics equations. This could be done by simply running this file as "main".
"""

if __name__ == "__main__":
     
    # -------------------------------------------------------------------------    
    # Define here the velocity vector. This is constant over time 
    # Wind tunnel approach
    # -------------------------------------------------------------------------
    
    u = 15.      # Component of free stream velocity
    v = None    # Component of lateral velocity
    w = -0.   # Component of vertical velocity
    
    # -------------------------------------------------------------------------
    # Define here the time array. Forces will be looped over each value
    # of time array.
    # -------------------------------------------------------------------------
    
    number_period = 1 # integer
    initial_time = 0
    final_time = number_period*(1/settings.frequency)
    time_steps = (50)*number_period
    time_array = np.linspace(initial_time, final_time, time_steps)
    
    # -------------------------------------------------------------------------
    # Call the FlappingForces function to run the Lifting Line code at 
    # This function is called for each time step
    # -------------------------------------------------------------------------
   
    lift_c = []
    drag = []
    moment = []
    tail = []
    moment_wing = []
    moment_tail = []
    moment_drag = []
    moment_lift = []
    for i in range (len(time_array)):
        
        [Fx, Fy, Fz, My, F_tail_tot, M_wing, M_tail, M_drag, M_lift] = FlappingForces(time_array[i], u, w, 
                                                                                        0,
                                                                                        0)
        lift_c.append(Fy)
        drag.append(Fz)
        moment.append(My)
        tail.append(F_tail_tot)
        moment_wing.append(M_wing)
        moment_tail.append(M_tail)
        moment_drag.append(M_drag)
        moment_lift.append(M_lift)
        print("Forces at time "+str(round(time_array[i], 3))+" ... done")
        
    Lift = np.asanyarray(lift_c)
    Drag = np.asanyarray(drag)
    Moment = np.asanyarray(moment)
    Tail = np.asanyarray(tail)
    Moment_wing = np.asanyarray(moment_wing)
    Moment_tail = np.asanyarray(moment_tail)
    Moment_drag = np.asanyarray(moment_drag)
    Moment_lift = np.asanyarray(moment_lift)
    
#    np.save('../utilities/forces_folder/lift_LL.npy', Lift)
#    np.save('../utilities/forces_folder/drag_LL.npy', Drag)

#    np.save('postprocessing/outfile_lift_u15', lift_c)

    import matplotlib.pyplot as plt
#
    fig = plt.figure()
    ax = fig.gca() 
    fig.suptitle('Lift comparison', fontsize=18)
    plt.xlabel('1/T', fontsize=14)
    plt.ylabel('L(t)', fontsize=14)
    ax.plot(time_array, Lift, '-', label = "Lift")
    ax.plot(time_array, Drag, '-', label = "Drag")
    ax.plot(time_array, Moment, '-', label = "Moment")

    ax.legend()
    ax.set_xlim(0, final_time)
    ax.grid(True)
#    
    fig1 = plt.figure()
    ax1 = fig1.gca() 
    fig1.suptitle('Moments', fontsize=18)
    plt.xlabel('1/T', fontsize=14)
    plt.ylabel('$M_i(t)$', fontsize=14)
    ax1.plot(time_array, Moment, '-', label = "Moment total")
    ax1.plot(time_array, Moment_tail, '-', label = "Moment tail")
    ax1.plot(time_array, Moment_wing, '-', label = "Moment wing")

    ax1.legend()
    ax1.set_xlim(0, final_time)
    ax1.grid(True)

    fig3 = plt.figure()
    ax3 = fig3.gca() 
    fig3.suptitle('Lift', fontsize=18)
    plt.xlabel('1/T', fontsize=14)
    plt.ylabel('$L(t)$', fontsize=14)
    ax3.plot(time_array, Lift, '-', label = "Lift")
    ax3.plot(time_array, Drag, '-', label = "Drag")

    ax3.legend()
    ax3.set_xlim(0, final_time)
    ax3.grid(True)

  
