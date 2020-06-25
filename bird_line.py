import numpy as np  # Import NumPy
import math as m
from wing_envelope import wing_envelope
from smooth_line import smooth_line
from smooth_quarterline import smooth_quarterline
from quarterchord_naive import quarterchord_naive
from improveline_iteration import improveline_iteration
import settings as settings

def bird_line(t, **kwargs):

    # Shoulder motion
    shoulder_x = bc.shoulder_x.motion_joint(t)
    shoulder_y = bc.shoulder_y.motion_joint(t)
    shoulder_z = bc.shoulder_z.motion_joint(t)
    # Elbow motion
    elbow_x = bc.elbow_x.motion_joint(t)
    elbow_y = bc.elbow_y.motion_joint(t)
    # Wrist motion
    wrist_y = bc.wrist_y.motion_joint(t)
    wrist_z = bc.wrist_z.motion_joint(t)

    """
    === Call of wing_envelope function. Given the kinematics, the wing shape is found ===
    """


    [leadingedge, trailingedge] =   wing_envelope(shoulder_x, shoulder_z, shoulder_y, elbow_y, elbow_x, wrist_y, wrist_z)

    nlete = 16

    leadingedge = smooth_line(leadingedge, nlete)
    trailingedge = smooth_line(trailingedge, nlete)

    [line, chord_leadingedge, chord_trailingedge] = quarterchord_naive(leadingedge, trailingedge)

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
        [line, chord_leadingedge, chord_trailingedge, a] = improveline_iteration(line,chord_leadingedge,
                                                        chord_trailingedge,leadingedge,trailingedge,it)

        chg = np.c_[chg, a]

        it = it + 1

    [lifting_line, chord_leadingedge, chord_trailingedge] = smooth_quarterline(line, chord_leadingedge, chord_trailingedge)

    """
    Output lifting_line, chord_leadingedge, chord_trailingedge CHECKED
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
    tArray = np.linspace(0, 1/f, 400)
    tip_line =[]
    root_line=[]
#    off_sweep = -np.deg2rad(26)
#    amp_sweep = np.deg2rad(20)
#    amplitude_shoulder = np.deg2rad(42)
#    amp_shoulder_x = 0.
#    settings.amplitude_shoulder_y = np.deg2rad(20)
    bc.shoulder_y = bc.Joint(-np.deg2rad(19), np.deg2rad(20), np.pi/2)
#    settings.offset_shoulder_y = -np.deg2rad(19)

    for i in range(len(tArray)):

        [lifting_line, updir, chord_direction,
         chord, line_left, up_direction_left,
         chord_direction_left, chord_left, dx] = bird_line(tArray[i])

#                                                             amp_shoulder_z=amplitude_shoulder,
#                                                             off_shoulder_y=off_sweep,
#                                                             amp_shoulder_y=amp_sweep)
        tip_line.append(lifting_line[:,-1])
        root_line.append(lifting_line[:,1])
    image_path = "/Users/gducci/UCL/MyPresentations/Confirmation/Figures"
    tip_line = np.asarray(tip_line)
    root_line=np.asarray(root_line)
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    wingframe_position = 0.05
    tip_line[:,2] = -tip_line[:,2]

    #   plt.ylim((0.,0.6))
    fig1 = plt.figure()
    ax1 = fig1.gca()
    plt.xlabel('$x_{airfoil}$', fontsize=18)
    plt.ylabel('$y_{airfoil}$', fontsize=18)
    plt.tick_params(labelsize=16)
    ax1.set_xlim(0.1, -0.3)  # decreasing time
#    ax1.set_xlim(0.45, -0.45)  # decreasing time


#    plt.scatter(-0., -0.0,color='black')
#    plt.scatter(wingframe_position, +0., color='black')
    switch_value = 0.
    supper = np.ma.masked_where(wingframe_position+tip_line[:,2] < switch_value, wingframe_position+tip_line[:,2])
    slower = np.ma.masked_where(wingframe_position+tip_line[:,2] > switch_value, wingframe_position+tip_line[:,2])
    smiddle = np.ma.masked_where((wingframe_position+tip_line[:,2] < switch_value) | (wingframe_position+tip_line[:,2] > switch_value), wingframe_position+tip_line[:,2])
    ax1.plot(supper, +0.+tip_line[:,1], slower, +0.+tip_line[:,1], smiddle, +0.+tip_line[:,1])
    plt.gcf().subplots_adjust(left=0.15)
    plt.gcf().subplots_adjust(bottom=0.15)

#    ax1.plot(wingframe_position+root_line[:,2],  root_line[:,1])
    plt.show()
#    ax1.grid()
    plt.savefig(image_path+'/Tip_trajectory.pdf', format = 'pdf')
#
#    function = off_sweep + amp_sweep*np.sin(((2*np.pi*f*tArray) + np.pi/2))
#    fig2 = plt.figure()
#    ax2 = fig2.gca()
#    plt.xlabel('$t$', fontsize=14)
#    plt.ylabel('$f(t) \ [deg]$', fontsize=14)
#    ax2.plot(tArray, np.rad2deg(function))
#    ax2.set_xlim(0,0.25)
#    ax2.plot(tArray, [np.mean(np.rad2deg(function))]*len(tArray), "--", color = "red", label = "Sweep offset")
#    plt.legend()
#    plt.savefig(image_path+'/Sweep_Function_0.eps', format = 'eps')

    fig3 = plt.figure(3)
    ax3 = fig3.gca(projection='3d')

    # Set axis label
    ax3.set_xlabel('$x$')
    ax3.set_ylabel('$y$')
    ax3.set_zlabel('$z$')


    # Plot the periodic orbit:
    ax3.plot(tip_line[:, 0],
            tip_line[:, 1],
            tip_line[:, 2],color = 'b')

