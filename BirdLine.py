import numpy as np  # Import NumPy
import math as m
from WingEnvelope import WingEnvelope 
from SmoothLine import SmoothLine
from SmoothQuarterLine import SmoothQuarterLine
from QuartChordNaive import QuartChordNaive
from ImproveLineIteration import ImproveLineIteration
import settings as settings
from KinematicsFunction import WingKinematics

def BirdLine(t, **kwargs):
# =============================================================================
#     Wrist joint default kinematics
# =============================================================================
    amp_shoulder_x=kwargs.get('amp_shoulder_x', 0.014)
    off_shoulder_x=kwargs.get('off_shoulder_x',0.2)
    phase_shoulder_x=kwargs.get('phase_shoulder_x',-np.pi/2)

    amp_shoulder_y=kwargs.get('amp_shoulder_y', np.pi/12)
    off_shoulder_y=kwargs.get('off_shoulder_y', -0.2 -np.pi/12)
    phase_shoulder_y=kwargs.get('phase_shoulder_y', np.pi/2)

    amp_shoulder_z=kwargs.get('amp_shoulder_z', np.deg2rad(45))
    off_shoulder_z=kwargs.get('off_shoulder_z',0.)
    phase_shoulder_z=kwargs.get('phase_shoulder_z', np.pi)
    
# =============================================================================
#     Elbow joint default kinematics
# =============================================================================
    off_elbow_y=kwargs.get('off_elbow_y', np.pi/6)
    amp_elbow_y=kwargs.get('amp_elbow_y', np.pi/6)
    phase_elbow_y=kwargs.get('phase_elbow_y', -np.pi/2)
    
    off_elbow_x=kwargs.get('off_elbow_x', 0.)
    amp_elbow_x=kwargs.get('amp_elbow_x', np.pi/6)
    phase_elbow_x=kwargs.get('phase_elbow_x', -np.pi/2)
    
# =============================================================================
#   Wrist default kinematics
# =============================================================================
    off_wrist_y = kwargs.get('off_wrist_y',-np.pi/6)
    amp_wrist_y = kwargs.get('amp_wrist_y',np.pi/6)
    phase_wrist_y = kwargs.get('phase_wrist_y',np.pi/2)

    off_wrist_z = kwargs.get('off_wrist_z',np.pi/12)
    amp_wrist_z = kwargs.get('amp_wrist_z',np.pi/12)
    phase_wrist_z = kwargs.get('phase_wrist_z',0.)
    
    shoulder_x = WingKinematics(off_shoulder_x, amp_shoulder_x, t, phase_shoulder_x)        # 0.014
    shoulder_y = WingKinematics(off_shoulder_y , amp_shoulder_y, t, phase_shoulder_y)              # -0.2 -np.pi/12, np.pi/12, t, np.pi/2
    shoulder_z = WingKinematics(off_shoulder_z, amp_shoulder_z, t, phase_shoulder_z)           # np.pi/4

    elbow_y = WingKinematics(off_elbow_y, amp_elbow_y, t, phase_elbow_y)
    elbow_x = WingKinematics(off_elbow_x, amp_elbow_x, t, phase_elbow_x)
    
    wrist_y = WingKinematics(off_wrist_y, amp_wrist_y, t, phase_wrist_y)
    wrist_z = WingKinematics(off_wrist_z, amp_wrist_z, t, phase_wrist_z)
   
#    shoulder_x = WingKinematics('shoulder', 'x', 0., 0., t,0.).motion       # 0,0,0
#    shoulder_y = WingKinematics('shoulder', 'y', -np.deg2rad(40), np.deg2rad(30), t, np.pi/2).motion    # -40 , 30, pi/2
#    shoulder_z = WingKinematics('shoulder', 'z', -np.deg2rad(0), np.deg2rad(60), t, np.pi).motion    # 0, 60
#
#    elbow_y = WingKinematics('elbow', 'y', np.deg2rad(45), np.deg2rad(45), t, -np.pi/2).motion  # 45, 45, -pi/2
#    elbow_x = WingKinematics('elbow', 'x', 0.,np.deg2rad(5), t, -np.pi/2).motion
#    
#    wrist_y = WingKinematics('wrist', 'y', -np.deg2rad(45), np.deg2rad(45), t, np.pi/2).motion # -45, 45, pi/2
#    wrist_z = WingKinematics('wrist', 'x', 0., 0., t, 0.).motion

    
    """
    === Call of WingEnvelope function. Given the kinematics, the wing shape is found ===
    """   

    
    [leadingedge, trailingedge] =   WingEnvelope(shoulder_x, shoulder_z, shoulder_y, elbow_y, elbow_x, wrist_y, wrist_z) 
    
    nlete = 16
        
    leadingedge = SmoothLine(leadingedge, nlete)
    trailingedge = SmoothLine(trailingedge, nlete)
    
    [line, chord_leadingedge, chord_trailingedge] = QuartChordNaive(leadingedge, trailingedge)
    
    """"
    ============= PLOT ROUTINE =============
    Here in the original code there is a function that prints out on the screen wing plots.
    It is now omitted, and there will be implemented later in a second time
    ========================================
    """
    
    tol = 0.1
    nmax = 10
    a = tol + 1
    chg = tol + 1
    it = 1
    
    while a > tol and it < nmax:
        [line, chord_leadingedge, chord_trailingedge, a] = ImproveLineIteration(line,chord_leadingedge,
                                                        chord_trailingedge,leadingedge,trailingedge,it)
        
        chg = np.c_[chg, a]
    
        it = it + 1
    
    [lifting_line, chord_leadingedge, chord_trailingedge] = SmoothQuarterLine(line, chord_leadingedge, chord_trailingedge)
    
    """
    Output lifting_line, chord_leadingedge, chord_trailingedge CHECKED. 
    """
    line_dummy = np.copy(lifting_line)
    nl = np.size(line[0])
    chord_direction = np.zeros((3,nl))
    updir = np.zeros((3,nl))
    chord = np.zeros((nl))
    
# =============================================================================
#     Evaluating the chord and its direction
#     For every slice, it's calculated the distance btw Lead. edge and Trail. 
#     edge (chord_distance), and then the modulus, so that the versor is 
#     identified.
# =============================================================================
    
    for i in range(nl):
        
        chord_distance = (chord_trailingedge[:,i] - chord_leadingedge[:,i]) + 1e-20
        chord[i] = np.linalg.norm(chord_distance)
        #chord = chord[:,np.newaxis]
        chord_direction[:,i] = chord_distance/chord[i]
        
        if i == 0:
            linevec = line[:,1] - line[:,0]
        elif i == (nl - 1):
            linevec = line[:,-1] - line[:,-2]
        else:
            linevec = line[:,i+1] - line[:,i-1]
            
        linevec = linevec/np.linalg.norm(linevec)
        updir[:,i] = np.cross(chord_direction[:,i], linevec)  # Different value in the last iteration
    
    updir_dummy = np.copy(updir)
    chord_direction_dummy = np.copy(chord_direction)
    chord_dummy = np.copy(chord)
    # Left Wing
    
    line_left = np.fliplr(line_dummy)
    up_direction_left = np.fliplr(updir_dummy)
    chord_direction_left = np.fliplr(chord_direction_dummy)
    chord_left = chord_dummy[::-1]
    
    line_left[0,:] = np.negative(line_left[0,:])
    chord_direction_left[0,:] = np.negative(chord_direction_left[0,:])
    up_direction_left[0,:] = np.negative(up_direction_left[0,:])
    
    sumvector = np.zeros((1,nl))
    for i in range(nl):
        if i < nl-1:
            sumvector[0,i] = np.sum((lifting_line[:,i+1] - lifting_line[:,i])**2)
        else:
            sumvector[0,i] = np.sum((lifting_line[:,-1] - lifting_line[:, -2])**2)
    
    
    dx = np.mean(np.sqrt(sumvector))
    
    return lifting_line, updir, chord_direction, chord, line_left, up_direction_left, chord_direction_left, chord_left, dx

if __name__ == "__main__":
    f = settings.frequency 
    tArray = np.linspace(0, 1/f, 200) 
    tip_line =[]
    off_sweep = -np.deg2rad(20)
    amp_sweep = np.deg2rad(40)
    amplitude_shoulder = np.deg2rad(46)
    for i in range(len(tArray)):

        [lifting_line, updir, chord_direction, chord, line_left, up_direction_left, chord_direction_left, chord_left, dx] = BirdLine(tArray[i],  amp_shoulder_z=amplitude_shoulder, 
                                                                                                                off_shoulder_y=off_sweep,
                                                                                                                amp_shoulder_y=amp_sweep)
        tip_line.append(lifting_line[:,-1])
           
    image_path = "/Users/gducci/UCL/MyPresentations/MeetingRenaud_1902_20"
    tip_line = np.asarray(tip_line)
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #   plt.ylim((0.,0.6))
    fig1 = plt.figure()
    ax1 = fig1.gca() 
    fig1.suptitle('Wing tip trajectory', fontsize=18)
    plt.xlabel('$z_{a}$', fontsize=14)
    plt.ylabel('$y_{a}$', fontsize=14)
    plt.scatter(-0., -0.0)
    plt.scatter(-0.04+tip_line[0,2], +0.+tip_line[0,1], color='green')
    switch_value = 0.
    
    supper = np.ma.masked_where(-0.04+tip_line[:,2] < switch_value, -0.04+tip_line[:,2])
    slower = np.ma.masked_where(-0.04+tip_line[:,2] > switch_value, -0.04+tip_line[:,2])
    smiddle = np.ma.masked_where((-0.04+tip_line[:,2] < switch_value) | (-0.04+tip_line[:,2] > switch_value), -0.04+tip_line[:,2])
    ax1.plot(supper, +0.+tip_line[:,1], slower, +0.+tip_line[:,1], smiddle, +0.+tip_line[:,1])
    plt.show()
#    plt.savefig(image_path+'/Tip_trajectory_0.eps', format = 'eps')
    
    function = off_sweep + amp_sweep*np.sin(((2*np.pi*f*tArray) + np.pi/2))
    fig2 = plt.figure()
    ax2 = fig2.gca()
    plt.xlabel('$t$', fontsize=14)
    plt.ylabel('$f(t) \ [deg]$', fontsize=14)
    ax2.plot(tArray, np.rad2deg(function))
    ax2.set_xlim(0,0.25)
    ax2.plot(tArray, [np.mean(np.rad2deg(function))]*len(tArray), "--", color = "red", label = "Sweep offset")
    plt.legend()
#    plt.savefig(image_path+'/Sweep_Function_0.eps', format = 'eps')

    fig3 = plt.figure()
    ax3 = fig3.gca()
    plt.xlabel('$t$', fontsize=14)
    plt.ylabel('$z(t)$', fontsize=14)
    ax3.plot(tArray, -0.04+tip_line[:,1])
    ax3.set_xlim(0,0.25)
    plt.legend()

#    for i in range(len(tip_line)):
#        if (-0.04+tip_line[i,2] < 0):
#            ax1.plot(-0.04+tip_line[i,2], +0.+tip_line[i,1], '-', color = "red")
#        if (-0.04+tip_line[i,2] > 0):
#            ax1.plot(-0.04+tip_line[i,2], +0.+tip_line[i,1], '-', color = "blue")
#    ax1.grid(True)



