import numpy as np  # Import NumPy
from BirdLine import BirdLine
from MergeLines import MergeLines
from Circulation import Circulation
import settings as settings

def FlappingForces(t, u, w, theta, q):
    cl_alpha = 2*np.pi
    rho = 1.225 # Air density
    
    U = np.array([0, w, u])
    #angular_velocity = np.array([q, 0, 0])
    dt = 1e-6
    t2 = t - dt
    wingframe_position = np.array([0, 0.1, -0.])
    wingframe_position_tail = np.array([0, 0.1, .4])
    
    # -------------------------------------------------------------------------
    # Initialization of the vector force. It must have the same dimension of the time vector
    # -------------------------------------------------------------------------
    Fx = 0
    Fy = 0
    Fz = 0
    My = 0
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
     line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(t)
    
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
     line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(t2)
    
    # -------------------------------------------------------------------------
    #  Geometry at time t2 = t-dt
    # -------------------------------------------------------------------------
    
    [line2, _, _, _] = MergeLines(line_right, up_direction_right,
    chord_direction_right, chord_right, line_left,up_direction_left, chord_direction_left, chord_left,nmid,x_ep)
    

    line_c = (line[:,0:-1] + line[:,1:])/2
    line_c2 = (line2[:,0:-1] + line2[:,1:])/2
    updir = (updir[:,0:-1]+updir[:,1:])/2
    chordir = (chordir[:,0:-1]+chordir[:,1:])/2
    chord = (chord[0:-1]+chord[1:])/2
    #velocity_q_induced = np.cross(angular_velocity, (line[2,:] + 0.1))

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
        
        [gamma_new, v_downwash] = Circulation(line[0,:], line[1,:], chord, angle_of_attack, cl_alpha, norm_velocity, gamma)
        
        error = np.mean(np.abs(gamma - gamma_new))/np.mean(np.abs(gamma_new))
        
        if error > error_p and count > 0:
            print("Error increasing after iteration ",count)
            break
        
        error_p = error
        
        gamma = gamma_new           
        
    gamma = np.reshape(gamma, np.size(gamma))
#    gamma[0:15] = np.flip(gamma[15::])
    
    local_velocity = velocity - v_downwash
    lever_arm = np.zeros_like(line_c)
    M = np.zeros_like(gamma)
    chord_tail = 0.15
    span_tail = 0.4
    Tail_Lift = 0.5*rho*(np.linalg.norm(U)**2)*chord_tail*cl_alpha*(np.arctan2(w, u))
    Tail_ds = np.linspace(0,span_tail, np.size(gamma))
    for j in range(np.size(gamma)):
        F = rho*gamma[j]*np.cross(local_velocity[:,j] , line_direction[:,j])
        lever_arm[:, j] = line_c[:, j] + wingframe_position
        F_tail_j = np.array([0, Tail_Lift*(Tail_ds[1]-Tail_ds[0]), 0.])
        F_moment = np.array([F[0]*ds[j], F[1]*ds[j], F[2]*ds[j]])
#        F_moment = np.array([F[0]*ds[j], F[1]*ds[j], F[2]*ds[j] + 0.5*((F[2]*ds[j]))])
#        F_moment_tail = np.array([0., 0.15*F[1]*ds[j], 0.])


        Fx  =   Fx + F[0]*ds[j]         # Component along wing span
        Fy  =   Fy + F[1]*ds[j]         # Component of Lift (perpendicular to the free stream)
        Fz  =   Fz + F[2]*ds[j]         # Component of drag (along free stream)
        M =  np.cross(lever_arm[:,j], F_moment) +  np.cross(wingframe_position_tail, F_tail_j)
        My = My + M[0]
        F_tail_tot = F_tail_tot + F_tail_j[1]
        
    return Fx, Fy, Fz, My, F_tail_tot



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
    w = 0.      # Component of vertical velocity
    
    # -------------------------------------------------------------------------
    # Define here the time array. Forces will be looped over each value
    # of time array.
    # -------------------------------------------------------------------------
    
    number_period = 1 # integer
    initial_time = 0
    final_time = number_period*(1/settings.frequency)
    time_steps = (100)*number_period
    time_array = np.linspace(initial_time, final_time, time_steps)
    
    # -------------------------------------------------------------------------
    # Call the FlappingForces function to run the Lifting Line code at 
    # This function is called for each time step
    # -------------------------------------------------------------------------
   
    lift_c = []
    drag = []
    moment = []
    tail = []
    for i in range (len(time_array)):
        
        [Fx, Fy, Fz, M, F_tail_tot] = FlappingForces(time_array[i], u, w, 0, 0)
        lift_c.append(Fy)
        drag.append(Fz)
        moment.append(M)
        print("Forces at time "+str(round(time_array[i], 3))+" ... done")
        
    Lift = np.asanyarray(lift_c)
    Drag = np.asanyarray(drag)
    Moment = np.asanyarray(moment)
    tail = np.asanyarray(tail)
#    np.save('../utilities/forces_folder/lift_LL.npy', Lift)
#    np.save('../utilities/forces_folder/drag_LL.npy', Drag)

#    np.save('postprocessing/outfile_lift_u15', lift_c)

    import matplotlib.pyplot as plt

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
