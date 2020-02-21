import numpy as np  # Import NumPy
from BirdLine import BirdLine
from MergeLines import MergeLines
from Circulation import Circulation
import settings as settings

def FlappingForces(t, u, w, q, theta):
    
    cl_alpha = 2*np.pi
    rho = 1.225 # Air density
    
    U = np.array([0, w, u])
    ang_velocity = np.array([q, 0, 0])

    dt = 1e-8
    t2 = t - dt
    
    
    """
    Initialization of the vector force. It must have the same dimension of the time vector
    """
    timelength = np.size(t)
    Fx = np.zeros(timelength)
    Fy = np.zeros(timelength)
    Fz = np.zeros(timelength)
    
    """
    For each time step t, the bird line is calculated, and therefore gamma and forces
    """
    x_ep = 0.05
    for i in range(timelength):
        
        #np.seterr(divide='ignore', invalid='ignore')

        
        """
        Evaluation of the geometry at time t
        """
        
        # -------------------------------------------------------------------------
        # Call bird line function, to extract the Lifting line of right and left wing
        # -------------------------------------------------------------------------
        
        [line_right, up_direction_right, chord_direction_right, chord_right, 
         line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(t[i])
        
        nmid = np.ceil(2*x_ep/dx)
        
        # -------------------------------------------------------------------------
        # Geometry at time t
        # -------------------------------------------------------------------------

        [line, chordir, updir, chord] = MergeLines(line_right, up_direction_right, chord_direction_right, chord_right, 
        line_left,up_direction_left, chord_direction_left, chord_left,nmid,x_ep)
        
        
        """
        Evaluation of the geometry at time t2 = t - dt
        """
        
        # -------------------------------------------------------------------------
        # Call bird line function, to extract the Lifting line of right and left wing
        # -------------------------------------------------------------------------
        
        [line_right, up_direction_right, chord_direction_right, chord_right, 
         line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(t2[i])

        # -------------------------------------------------------------------------
        # Geometry at time t2 = t - dt
        # -------------------------------------------------------------------------

        [line2, _, _, _] = MergeLines(line_right, up_direction_right, chord_direction_right, chord_right, 
        line_left,up_direction_left, chord_direction_left, chord_left,nmid,x_ep)
        
        print("Line t-step ", i, "OK")
        line_c = (line[:,0:-1] + line[:,1:])/2          # Mean of discretised point at time t
        line_c2 = (line2[:,0:-1] + line2[:,1:])/2       # Mean of discretised point at time t2 = t-dt
        updir = (updir[:,0:-1]+updir[:,1:])/2
        chordir = (chordir[:,0:-1]+chordir[:,1:])/2
        chord = (chord[0:-1]+chord[1:])/2
        # -------------------------------------------------------------------------
        # u_inf = free stream velocity, equal over each profile (points)
        # -------------------------------------------------------------------------
        #angular_velocity = np.array([ang_velocity for j in range(np.size(line_c[0]) + 1)]).T
        u_inf = np.array([U for j in range(np.size(line_c[0]))]).T
        #velocity_q_induced[:,:] = np.cross(angular_velocity[:,:], line[:,:])

        # -------------------------------------------------------------------------
        # ???? Double check velocity, because the flapping velocity has a direction!
        # -------------------------------------------------------------------------

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
            [gamma_new, vn] = Circulation(line[0,:], line[1,:], chord, angle_of_attack, cl_alpha, np.linalg.norm(U), gamma)
            
            error = np.mean(np.abs(gamma - gamma_new))/np.mean(np.abs(gamma_new))

            if error > error_p and count > 0:
                print("Error increasing for t-step ", i, "after iteration", count)
                break
            
            error_p = error
            
            gamma = gamma_new
            
            v_downwash = vn
            
        local_velocity = velocity - v_downwash
        
        for j in range(np.size(gamma)):
            F = rho*gamma[0,j]*np.cross(local_velocity[:,j] , line_direction[:,j])
            Fx[i]  =   Fx[i] + F[0]*chord[j]*ds[j]       # Component along wing span
            Fy[i]  =   Fy[i] + F[1]*chord[j]*ds[j]       # Component of Lift (perpendicular to the free stream)
            Fz[i]  =   Fz[i] + F[2]*chord[j]*ds[j]       # Component of drag (along free stream)
        
    
#    from pylab import subplot, plot, xlabel, ylabel, show, scatter
#    import matplotlib.pyplot as plt
#    #
#    plt.scatter(t, gamma)
#    ylabel('x(t)')
#    
#    xlabel('t (s)')
#    ylabel('v(t)')
#    
#    
#    plt.show()    
    
    return Fx, Fy, Fz


if __name__ == "__main__":
    
   f = settings.frequency
   
   period = 1/f
   
   t_disc = np.linspace(0,period,50)
   
   [Fx_0, Fy_0, Fz_0] = FlappingForces(t_disc, 10, 0, 0, 0)
   print(Fy_0)
   t = len(Fx_0)

   import matplotlib.pyplot as plt
   fig, ax = plt.subplots()
   ax.plot(t_disc, Fy_0)
    
   ax.set(xlabel='t/T', ylabel='$L_{nc}$',
           title='lift')
   ax.grid()

   
  