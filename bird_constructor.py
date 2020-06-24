import numpy as np
import math as m
from WingEnvelope import WingEnvelope
from SmoothLine import SmoothLine
from SmoothQuarterLine import SmoothQuarterLine
from QuartChordNaive import QuartChordNaive
from ImproveLineIteration import ImproveLineIteration
import settings as settings
from KinematicsFunction import WingKinematics

frequency  = 4
omega = frequency*2*np.pi

class Joint:
    omega = frequency*2*np.pi

    def __init__(self, offset, amplitude, phase):
        self.offset = offset
        self.amplitude = amplitude
        self.phase = phase


    def motion_joint(self, t):

        motion = self.offset + self.amplitude*np.sin(omega*t + self.phase)
        return motion

class Shoulder(Joint):

    def __init__(self, offset, amplitude, phase, axis):

        super().__init__(offset, amplitude, phase)
        self.axis = axis

class Elbow(Joint):

    def __init__(self, offset, amplitude, phase, axis):

        super().__init__(offset, amplitude, phase)
        self.axis = axis

class Wrist(Joint):

    def __init__(self, offset, amplitude, phase, axis):

        super().__init__(offset, amplitude, phase)
        self.axis = axis

class Bird:

    def __init__(self, offset, amplitude, phase):
        self.shoulder_x = Shoulder(offset, amplitude, phase, 'x')
        self.shoulder_y = Shoulder(offset, amplitude, phase, 'y')
        self.shoulder_z = Shoulder(offset, amplitude, phase, 'z')
        self.elbow_x =  Elbow(offset, amplitude, phase, 'x')
        self.elbow_y = Elbow(offset, amplitude, phase, 'y')
        self.wrist_y = Wrist(offset, amplitude, phase, 'y')
        self.wrist_z = Wrist(offset, amplitude, phase, 'z')

    def BirdLine(self, t):

        # Shoulder motion
        shoulder_x = self.shoulder_x.motion_joint(t)
        shoulder_y = self.shoulder_y.motion_joint(t)
        shoulder_z = self.shoulder_z.motion_joint(t)
        # Elbow motion
        elbow_x = self.elbow_x.motion_joint(t)
        elbow_y = self.elbow_y.motion_joint(t)
        # Wrist motion
        wrist_y = self.wrist_y.motion_joint(t)
        wrist_z = self.wrist_z.motion_joint(t)
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
            updir[:,i] = np.cross(chord_direction[:,i], linevec)

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

# Initialization of the kinematics
#shoulder_x = Joint(0.2, 0.014, -np.pi/2)
shoulder_x = Shoulder(0.2, 0.014, -np.pi/2, 'x')
shoulder_y = Shoulder(-np.deg2rad(19), np.deg2rad(20), np.pi/2, 'y')
shoulder_z = Shoulder(0., np.deg2rad(42), np.pi, 'z')

elbow_x = Elbow(0., np.pi/6, -np.pi/2, 'x')
elbow_y = Elbow(np.pi/6, np.pi/6, -np.pi/2, 'y')

wrist_y = Wrist(-np.pi/6, np.pi/6, np.pi/2, 'y')
wrist_z = Wrist(0., 0., 0., 'z')
