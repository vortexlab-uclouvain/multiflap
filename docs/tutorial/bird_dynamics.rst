Bird model dynamics
===================

This tutorial reproduces step-by-step the results obtained analysing the bird dynamics.  

Imposing the bird kinematics
****************************

The bird wing is a poliarticulated rigid body, with three joints: shoulder, elbow and wrist. Each joint motion is governed by a periodic function of the form:

.. math::
   :label: joint_kinematics

   q_i(t) = q_{0,i} + A_i \sin\left(\omega t + \phi_{0,i}\right)
  
where :math:`q_{i}` is the rotation of the joint with respect to the generic axis :math:`i`.

The module ``multiflap/aero_package/kinematics_constructor.py`` allows to impose the kinamatics on each wing joint.


.. toggle-header::
    :header: ``kinematics_constructor.py`` **Show code**

            .. code-block:: python

                import numpy as np
                import math as m

                frequency  = 4
                omega = frequency*2*np.pi

                class Joint:
                    omega = frequency*2*np.pi

                    def __init__(self, offset=None, amplitude=None, phase=None):
                        self.offset = offset
                        self.amplitude = amplitude
                        self.phase = phase

                    def motion_joint(self, t):

                        motion = self.offset + self.amplitude*np.sin(omega*t + self.phase)
                        return motion

                class Shoulder():
                    def __init__(self, axis_x=None, axis_y=None, axis_z=None):
                        self.axis_x  = Joint()
                        self.axis_y  = Joint()
                        self.axis_z  = Joint()
                        if isinstance(axis_x, Joint):
                            self.axis_x = axis_x
                        if isinstance(axis_y, Joint):
                            self.axis_y = axis_y
                        if isinstance(axis_z, Joint):
                            self.axis_z = axis_z

                class Elbow():
                    def __init__(self, axis_x=None, axis_y=None):
                        self.axis_x  = Joint()
                        self.axis_y  = Joint()
                        if isinstance(axis_x, Joint):
                            self.axis_x = axis_x
                        if isinstance(axis_y, Joint):
                            self.axis_y = axis_y

                class Wrist():
                    def __init__(self, axis_y=None, axis_z=None):
                        self.axis_y  = Joint()
                        self.axis_z  = Joint()
                        if isinstance(axis_y, Joint):
                            self.axis_y = axis_y
                        if isinstance(axis_z, Joint):
                            self.axis_z = axis_z

The bird object
***************

Once the kinematics is imposes, the object **bird** can be crated. This object is the input for the aerodynamic solver because it takes the kinematics of each joint, and exctracts the wing envelope and subsequently the lifting line at each time step.

The bird object takes as arguments the kinematics previously generated. The module is reported below.

.. toggle-header::
    :header: ``bird_model.py`` **Show code**

            .. code-block:: python


                import numpy as np
                from .kinematics_constructor import Joint, Shoulder, Elbow, Wrist
                from .line_functions import smooth_line, smooth_quarterline, quarterchord_naive, improveline_iteration
                from .RotationMatrix import RotationMatrix
                from .CompMatrix import CompMatrix
                from .Plumes import plumes
                from collections import namedtuple
                from odes.bird_dynamics import dynamics, get_aeroforces, get_stability_matrix


                class BirdModel:

                    def __init__(self, shoulder=None, elbow=None, wrist=None):

                        #Parameters:
                        self.dimensions = 4
                        self.g = 9.81
                        self.mass = 1.2
                        self.frequency = 4
                        self.wingframe_position = np.array([0, -0.0, -0.05])
                        self.wingframe_position_tail = np.array([0, -0., .3])
                        self.tail_length = 0.25
                        self.tail_chord = 0.15
                        self.tail_opening = 0.
                        self.shoulder = Shoulder()
                        self.wrist = Wrist()
                        self.elbow = Elbow()

                        if isinstance(shoulder, Shoulder):
                            self.shoulder = shoulder
                        if isinstance(wrist, Wrist):
                            self.wrist = wrist
                        if isinstance(elbow, Elbow):
                            self.elbow = elbow


                    def get_wingstate(self, t):

                        state_shoulder_x = self.shoulder.axis_x.motion_joint(t)
                        state_shoulder_y = self.shoulder.axis_y.motion_joint(t)
                        state_shoulder_z = self.shoulder.axis_z.motion_joint(t)
                        # Elbow motion
                        state_elbow_x = self.elbow.axis_x.motion_joint(t)
                        state_elbow_y = self.elbow.axis_y.motion_joint(t)
                        # Wrist motion
                        state_wrist_y = self.wrist.axis_y.motion_joint(t)
                        state_wrist_z = self.wrist.axis_z.motion_joint(t)

                        wing_state = [state_shoulder_x, state_shoulder_z, state_shoulder_y, state_elbow_y, state_elbow_x, state_wrist_y, state_wrist_z]
                        #return ws
                        return wing_state

                    def get_wingenvelope(self, wing_state):
                        shoulder_x, shoulder_z, shoulder_y, elbow_y, elbow_x, wrist_y, wrist_z = wing_state
                        [plumesrot, plumesvec] = plumes()
                        dim_plumesrot0 = np.size(plumesrot[0])
                        """
                        Assignment of the dimensions of the skeleton and the frame position
                        """

                        # Versor along y axis
                        ey = [0, 1, 0]

                        # Lenght of the first part of the wing skeleton
                        l1 = .134

                        # Lenght of the second part of the wing skeleton
                        l2 = .162

                        # Lenght of the second third of the wing skeleton (l3 is used to connect l1 and l2)
                        l4 = .084

                        # Origin of the reference frame
                        origin = [0, 0, 0]

                        # Respective position of the end of the three parts of the skeleton
                        end1 = np.array([l1, 0, 0])
                        end2 = np.array([l2, 0, 0])
                        end4 = np.array([l4, 0, 0])

                        vec30 = np.array(end1) + np.array(end2) - np.array(origin) # Probably unused

                        """
                        Assignment of the points over the three parts of the skeleton, arm, forearm, hand respoectively
                        """
                        # Number of points over the arm
                        npoints_arm = 1

                        # Number of points over the forearm
                        npoints_forearm = 1

                        # Number of points over the hand
                        npoints_hand = 2

                        # Total number of point
                        nslice = npoints_arm + npoints_forearm + npoints_hand

                        # Distances from point to point for each part of the skeleton
                        delta_arm = l1/npoints_arm
                        delta_forearm = l2/npoints_forearm
                        delta_hand = l4/(npoints_hand-1)


                        # Routine for the xslice vector (values checked with Matlab code). Maybe write a function!
                        xslice = np.zeros(nslice)
                        for i in range(npoints_arm):
                            xslice[i] = 0 + (i*delta_arm)
                        for i in range(npoints_forearm):
                            xslice[i + npoints_arm] = l1 + (i*delta_forearm)
                        for i in range(npoints_hand):
                            xslice[i + npoints_arm + npoints_forearm] = l1 + l2 + (i*delta_hand)

                        # Index of feathers starting point (ask it for a more clear comment). CHECKED!
                        index = np.zeros(7)
                        index[0] = np.size(xslice)
                        index[1] = np.nonzero(xslice>=l1+l2+l4/2)[0][0]
                        index[2] = np.nonzero(xslice>=l1+l2)[0][0]
                        index[3] = np.nonzero(xslice>=l1+l2/2)[0][0]
                        index[4] = np.nonzero(xslice>=l1)[0][0]
                        index[5] = np.nonzero(xslice>=l1/2)[0][0]
                        index[6] = 1
                        number_ellipses = 2
                        s = np.linspace(0,2*np.pi,number_ellipses + 1) 
                        s = s[0:number_ellipses]
                        skz = .01
                        sky = .01
                        z = skz*np.cos(s)
                        y = sky*np.sin(s)
                        points = np.zeros((3 , nslice*number_ellipses))

                        w1 = np.zeros(nslice*number_ellipses)
                        w2a = np.zeros(nslice*number_ellipses)
                        w2b = np.zeros(nslice*number_ellipses)
                        w3 = np.zeros(nslice*number_ellipses)
                        w4 = np.zeros(nslice*number_ellipses)
                        for j in range(nslice):
                            if xslice[j] < l1:
                                for k in range(number_ellipses):
                                    points[:,(j)*number_ellipses+k] = [xslice[j], y[k], z[k]]

                                facb = (z/skz + 1)/2
                                lrefcd = min(l1,l2)

                                if xslice[j] < l1-lrefcd/2:
                                    faccd1 = 1
                                else:
                                    faccd1 = 0.5*(1-np.sin(np.pi*(xslice[j]-l1)/lrefcd))

                                faccd2 = 1-faccd1
                                fac2a = 1
                                fac2b = 0
                                w1[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd1
                                w2a[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2a
                                w2b[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2b
                                w3[j*number_ellipses:(j+1)*number_ellipses] = (1-facb)
                                w4[j*number_ellipses:(j+1)*number_ellipses] = 0

                            elif xslice[j] < l1 + l2:
                                for k in range(number_ellipses):
                                    points[:,(j)*number_ellipses+k] = [xslice[j], y[k], z[k]]

                                facb = (z/skz + 1)/2
                                lrefcd = min(l1,l2)
                                if xslice[j] < l1+lrefcd/2:
                                    faccd2 = 0.5*(1+np.sin(np.pi*(xslice[j]-l1)/lrefcd))
                                else:
                                    faccd2 = 1

                                faccd1 = 1-faccd2
                                fac2b = (xslice[j]-l1)/l2
                                fac2a = 1-fac2b
                                lrefpg = min(l2,l4)

                                if xslice[j] < l1+l2-lrefpg/2:
                                    facpg2 = 1
                                else:
                                    facpg2 = 0.5*(1-np.sin(np.pi*(xslice[j]-(l1+l2))/lrefpg))

                                facpg4 = 1-facpg2

                                w1[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd1*facpg2
                                w2a[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2a*facpg2
                                w2b[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2b*facpg2
                                w3[j*number_ellipses:(j+1)*number_ellipses] = (1-facb)*facpg2
                                w4[j*number_ellipses:(j+1)*number_ellipses] = facpg4
                            else:
                                for k in range(number_ellipses):
                                    points[:,(j)*number_ellipses+k] = [xslice[j], y[k], z[k]]

                                facb = (z/skz + 1)/2
                                fac2a = 0
                                fac2b = 1
                                lrefpg = min(l2,l4)

                                if xslice[j] < l1+l2+lrefpg/2:
                                    facpg4 = 0.5*(1+np.sin(np.pi*(xslice[j]-(l1+l2))/lrefpg))
                                else:
                                    facpg4 = 1

                                facpg2 = 1-facpg4

                                w1[j*number_ellipses:(j+1)*number_ellipses] = 0
                                w2a[j*number_ellipses:(j+1)*number_ellipses] = facb*facpg2*fac2a
                                w2b[j*number_ellipses:(j+1)*number_ellipses] = facb*facpg2*fac2b
                                w3[j*number_ellipses:(j+1)*number_ellipses] = (1-facb)*facpg2
                                w4[j*number_ellipses:(j+1)*number_ellipses] = facpg4



                        npts = nslice*number_ellipses

                        pts = points

                        """
                        Rotation routine for the three parts of the skeleton (shoulder, elbow, hand)
                        """
                        rot1X = RotationMatrix(shoulder_x,'x')
                        rot1Z = RotationMatrix(shoulder_z,'z')
                        rot1Y = RotationMatrix(shoulder_y,'y')
                        rot1 = np.array(rot1X).dot(np.array(rot1Z)).dot(np.array(rot1Y))
                        ey1 = np.array(rot1).dot(np.array(ey))

                        rot2relY = RotationMatrix(elbow_y,'y')
                        rot2relX = RotationMatrix(elbow_x,'x')
                        rot2b = np.array(rot1).dot(np.array(rot2relY)).dot(np.array(rot2relX))
                        rot2a = np.array(rot1).dot(np.array(rot2relY))

                        rot4relY = RotationMatrix(wrist_y,'y')
                        rot4relZ = RotationMatrix(wrist_z,'z')
                        rot4 = np.array(rot2b).dot(np.array(rot4relY)).dot(np.array(rot4relZ))


                        """
                        Evaluation of the bones position
                        """
                        endd1 = np.array(rot1).dot(np.array(end1))
                        endd2 = np.array(rot2a).dot(np.array(end2))+np.array(endd1)
                        endd4 = np.array(rot4).dot(np.array(end4))+np.array(endd2)
                        vec3 = np.array(endd2)-np.array(origin)

                        ex3 = vec3/np.linalg.norm(vec3)
                        ey3 = ey1
                        ez3 = np.array([ex3[1]*ey3[2]-ex3[2]*ey3[1], ex3[2]*ey3[0]-ex3[0]*ey3[2], ex3[0]*ey3[1]-ex3[1]*ey3[0]])
                        rot3 = np.array([ex3,ey3,ez3])
                        rot3 = rot3.transpose()
                        kx = np.linalg.norm(vec3)/np.linalg.norm(vec30)
                        comp3 = CompMatrix(kx,'x')


                        """
                        Routine to find skin points
                        """

                        for j in range(npts):
                            pt = np.array(points[:,j])
                            pt1 = rot1.dot(pt)
                            pt2a = rot2a.dot((pt-end1)) + endd1
                            pt2b = rot2b.dot((pt-end1)) + endd1
                            pt3 = rot3.dot(comp3).dot(pt)
                            pt4 = rot4.dot((pt-(end1+end2))) + endd2

                            pts[:,j] = pt1*w1[j] + pt2a*w2a[j] + pt2b*w2b[j] + pt3*w3[j] + pt4*w4[j]

                        """
                        Routine to find main feathers
                        """
                        vecpss = np.zeros((3,dim_plumesrot0))
                        startp = np.zeros((3,np.size(plumesrot[0])))

                        for i in range(np.size(plumesrot[0])):
                            if i == 0:
                                prerot = rot4
                            elif i == 1:
                                prerot = rot4.dot(RotationMatrix(wrist_z/4,'x'))
                            elif i == 2:
                                prerot = rot2b.dot(RotationMatrix(wrist_y/2,'y')).dot(RotationMatrix(-wrist_z/2,'x'))
                            elif i == 3:
                                prerot = rot2a.dot(RotationMatrix(elbow_x/2,'x')).dot(RotationMatrix(wrist_y/2,'y')).dot(RotationMatrix(-wrist_z/4,'x'))
                            elif i == 4:
                                prerot = rot1.dot(RotationMatrix(elbow_y/2,'y'))
                            elif i == 5:
                                prerot = rot1.dot(RotationMatrix(elbow_y/4,'y'))
                            elif i == 6:
                                prerot = rot1

                            rY = RotationMatrix(plumesrot[0,i],'y')
                            rX = RotationMatrix(plumesrot[1,i],'x')
                            rot = prerot.dot(rY).dot(rX)
                            vecpss[:,i] = np.dot(rot,plumesvec[:,i])
                            element = index[i]
                            column = int((index[i]-1)*(number_ellipses))
                            startp[0,i] = pts[0][column]
                            startp[1,i] = pts[1][column]
                            startp[2,i] = pts[2][column]

                        """
                        Interpolation Routine
                        """
                        main1 = ((xslice[xslice >= l1+l2+l4/2]-(l1+l2+l4/2))/(l4/2))
                        main1 = main1[:,np.newaxis].T
                        main2 = ((xslice[(xslice >= l1+l2) & (xslice < l1+l2+l4/2)]-(l1+l2))/(l4/2))
                        main2 = main2[:,np.newaxis].T
                        abras1 = ((xslice[(xslice >=l1+l2/2) & (xslice <l1+l2)]-(l1+l2/2))/(l2/2))
                        abras1 = abras1[:,np.newaxis].T
                        abras2 = ((xslice[(xslice >=l1) & (xslice <l1+l2/2)]-(l1))/(l2/2))
                        abras2 = abras2[:,np.newaxis].T
                        bras1 = ((xslice[(xslice >=l1/2) & (xslice <l1)]-(l1/2))/(l1/2))
                        bras1 = bras1[:,np.newaxis].T
                        bras2 = (xslice[xslice <l1/2]/(l1/2))
                        bras2 = bras2[:,np.newaxis].T

                        """
                        Feathers Routine
                        """
                        pm1 = vecpss[:,0][:,np.newaxis].dot(main1)+vecpss[:,1][:,np.newaxis].dot(1-main1)
                        pm2 = vecpss[:,1][:,np.newaxis].dot(main2)+vecpss[:,2][:,np.newaxis].dot(1-main2)
                        pa1 = vecpss[:,2][:,np.newaxis].dot(abras1)+vecpss[:,3][:,np.newaxis].dot(1-abras1)
                        pa2 = vecpss[:,3][:,np.newaxis].dot(abras2)+vecpss[:,4][:,np.newaxis].dot(1-abras2)
                        pb1 = vecpss[:,4][:,np.newaxis].dot(bras1)+vecpss[:,5][:,np.newaxis].dot(1-bras1)
                        pb2 = vecpss[:,5][:,np.newaxis].dot(bras2)+vecpss[:,6][:,np.newaxis].dot(1-bras2)

                        """
                        Store the feathers in a single vector following the structure:

                                [pb2(0,0),..,pb2(0,n),pb1(0,0),..,pb1(0,n),pa2(0,0),..,pa2(0,k),pa1(0,0),..,pa1(0,k),pm2(0,0),..,pm2(0,j),pm1(0,0),..,pm1(0,j)]
                        vecps = [pb2(1,0),..,pb2(1,n),pb1(1,0),..,pb1(1,n),pa2(1,0),..,pa2(1,k),pa1(1,0),..,pa1(1,k),pm2(1,0),..,pm2(1,j),pm1(1,0),..,pm1(1,j)]
                                [pb2(2,0),..,pb2(2,n),pb1(2,0),..,pb1(2,n),pa2(2,0),..,pa2(2,k),pa1(2,0),..,pa1(2,k),pm2(2,0),..,pm2(2,j),pm1(2,0),..,pm1(2,j)]

                        Matrix concatenation
                        """
                        vecps = np.c_[pb2, pb1, pa2, pa1, pm2, pm1]


                        trailingedge = np.zeros((3,nslice))

                        xpt = np.zeros(3)
                        ypt = np.zeros(3)
                        zpt = np.zeros(3)

                        for j in range(nslice):
                            xpt = np.array([pts[0,(j)*number_ellipses], pts[0,j*number_ellipses+1], pts[0,(j)*number_ellipses]])
                            ypt = np.array([pts[1,(j)*number_ellipses], pts[1,j*number_ellipses+1], pts[1,(j)*number_ellipses]])
                            zpt = np.array([pts[2,(j)*number_ellipses], pts[2,j*number_ellipses+1], pts[2,(j)*number_ellipses]])

                            vecp = np.array(vecps[:,j])

                            trailingedge[0,j] = xpt[0]+vecp[0]
                            trailingedge[1,j] = ypt[0]+vecp[1]
                            trailingedge[2,j] = zpt[0]+vecp[2]

                        leadingedge = np.c_[pts[:,int(number_ellipses/2)::int(number_ellipses)],trailingedge[:,-1]]
                        return leadingedge, trailingedge
                    def get_liftingline(self, leadingedge, trailingedge):
                        """
                        === Call of wing_envelope function. Given the kinematics, the wing shape is found ===
                        """

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



                    def merge_lines(self, line_package):
                        line_r,updir_r,chordir_r,chord_r,line_l,updir_l,chordir_l,chord_l, dx = line_package
                        x_ep = 0.05
                        nmid = int(np.ceil(2*x_ep/dx))
                        line_r[0,:] = line_r[0,:] + x_ep
                        line_l[0,:] = line_l[0,:] - x_ep

                        xmid = np.linspace(line_l[0,-1],line_r[0,0],nmid)
                        length_xmid = np.size(xmid)
                        line_mid = np.array([xmid, line_r[1,0]*np.ones(np.size(xmid)), line_r[2,0]*np.ones(np.size(xmid))])
                        chordir_mid = np.array([chordir_r[:,0] for i in range(np.size(xmid))]).T
                        updir_mid = np.array([updir_r[:,0] for i in range(np.size(xmid))]).T
                        chord_mid = np.zeros(length_xmid)
                        chord_mid[:] = chord_r[0]*np.ones(length_xmid)

                        line = np.c_[line_l,line_mid[:,1:-1],line_r]
                        chordir = np.c_[chordir_l, chordir_mid[:,1:-1], chordir_r]
                        updir = np.c_[updir_l, updir_mid[:,1:-1], updir_r]
                        chord = np.concatenate([chord_l, chord_mid[1:-1]/1e4, chord_r])
                        return line, chordir, updir, chord

                    # calling methods of BirdModel that are defined in other files
                    get_aeroforces = get_aeroforces
                    dynamics = dynamics
                    get_stability_matrix = get_stability_matrix


.. Note::
   Kinematics and bird constructors have been briefly introduced because they will be called from the user in the ``main.py`` file.

We are now able to create the user files.

Input file
**********

As all the other input files, this has to respect the path ``multiflap/odes/bird_dynamics.py``. The equations of motion of the flier restricted on the longitudinal plane are:

.. math::
   :label: bird_odes

   \dot{u} &= -qw - g\sin \theta + \frac{1}{m_b}{F_{x'}(\mathbf{x}(t),  t)}\\ 
   \dot{w} &= qu + g\cos \theta +\frac{1}{m_b} F_{z'}(\mathbf{x}(t),  t)\\
   \dot{q} &=\frac{1}{I_{yy}}M_{y'}(\mathbf{x}(t),  t)\\
   \dot{\theta} &= q

and therefore the time-dependent stability matrix:

.. math::
   :label: stability_matrix_bird


   \mathbb{A}(\mathbf{x}(t), t) = 
   \begin{pmatrix}
   \frac{1}{m}\frac{\partial{F_{x'}}}{\partial{u}}  & \frac{1}{m}\frac{\partial{F_{x'}}}{\partial{w}} - q& - w + \frac{1}{m}\frac{\partial{F_{x'}}}{\partial{q}}& -g\cos{\theta} + \frac{1}{m}\frac{\partial{F_{x'}}}{\partial{\theta}}\\
   q + \frac{1}{m}\frac{\partial{F_{z'}}}{\partial{u}}  & \frac{1}{m}\frac{\partial{F_{z'}}}{\partial{w}}&  u + \frac{1}{m}\frac{\partial{F_{z'}}}{\partial{q}} & g\sin{\theta} + \frac{1}{m}\frac{\partial{F_{z'}}}{\partial{\theta}}\\
   \frac{1}{I_{yy}}\frac{\partial{M_{y'}}}{\partial{u}}  & \frac{1}{I_{yy}}\frac{\partial{M_{y'}}}{\partial{w}}& \frac{1}{I_{yy}}\frac{\partial{M_{y'}}}{\partial{q}}  & \frac{1}{I_{yy}}\frac{\partial{M_{y'}}}{\partial{\theta}}\\
   0 & 0 &1 & 0 \\
   \end{pmatrix}

The aerodynamic derivatives are computed at every time step from the ``get_stability_matrix`` function.

The input file is ``multiflap/odes/bird_dynamics.py``. This set of ODEs is coupled with the aerodynamics. To use the aerodynamic solver, the realative package ``aero_package`` needs to be loaded.

At every time step the position of the wing and the extraction of the lifting line is calculated, based on the bird object passed from the main file.


.. code-block:: python

        import numpy as np
        from aero_package.flapping_forces import get_flappingforces
        import scipy.integrate as ode
        import time
        from math import sin, cos
        dim = 4
        
        def dynamics(self, x0, t):
                """
                ODE system for the Equation of Motion on the longitudinal plane.
                This function will be passed to the numerical integrator

                Inputs:
                    x0: initial values [u, w, q, theta]
                    t: time

                Outputs:
                    x_dot: velocity vector
                """
                #Read inputs:
                u, w, q, theta  = x0  # Read state space points

                # get the aereodynamic forces at time t
                [_, Fy, Fz, My, F_tail, _, _, _, _] = self.get_aeroforces(x0, t)

                # bird body dynamics
                dudt = -q*w - self.g*np.sin(theta) - Fz/self.mass
                dwdt = q*u + self.g*np.cos(theta) - Fy/self.mass - F_tail/self.mass
                dqdt =  My/0.1
                dthetadt = q

                # collecting x_dot components in a single array:
                x_dot = np.array([dudt, dwdt, dqdt, dthetadt], float)

                return x_dot

        def get_aeroforces(self, x0, t):
                """
                Returns the aereodynamic forces at a time t
                1. Get the wing state (angles of each joint at time t)
                2. Based on the joint angles, exctracts the wing envelope
                3. From the envelope, exctracts the lifing line
                4. The lifting line is mirrored for the symmetric wing
                5. Lifting line properties and bird velocities are
                    passed to the liftin line solver

                Inputs:
                    x0: initial values [u, w, q, theta]
                    t: time

                Outputs:
                    aereodynamic_forces : components of the aereod. forces and moments
                """
                dt = 0.25*1e-4
                wing_state = self.get_wingstate(t)
                wing_state_dt = self.get_wingstate(t-dt)

                [leadingedge,
                 trailingedge] = self.get_wingenvelope(wing_state)

                [leadingedge_dt,
                 trailingedge_dt] = self.get_wingenvelope(wing_state_dt)

                lifting_line =  self.get_liftingline(leadingedge,
                                                     trailingedge)

                lifting_line_dt =  self.get_liftingline(leadingedge_dt,
                                                        trailingedge_dt)

                [line, chordir, updir, chord] = self.merge_lines(lifting_line)

                [line_dt, _, _, _] = self.merge_lines(lifting_line_dt)

                v_kinematics = get_kinematicvelocity(line, line_dt, dt)

                # call the aereodynamic solver for the lifting line get_flappingforces
                [Fx,
                 Fy,
                 Fz,
                 My,
                 F_tail,
                 M_wing,
                 M_tail,
                 M_drag,
                 M_lift] = get_flappingforces(x0, v_kinematics,
                                              line, chordir, updir, chord)

                aereodynamic_forces = [Fx, Fy, Fz, My,
                                       F_tail, M_wing, M_tail, M_drag, M_lift]
                return aereodynamic_forces


        def get_kinematicvelocity(line, line_dt, dt):

               line_c = line[:, 1:-1:2]
               line_c2 = line_dt[:, 1:-1:2]
               v_kin =  (line_c - line_c2)/dt

               return v_kin

        def get_stability_matrix(self, x0, t):

            """
            Stability matrix of the ODE system for the longitudinal plane

            Inputs:
                x0: initial condition [u_0, w_0, q_0, theta_0]
            Outputs:
                A: Stability matrix evaluated at x0. (dxd) dimension
                A[i, j] = dv[i]/dx[j]
            """
            m = self.mass
            g = self.g
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

            A = np.array([[-dFzu_dU/m, -q - dFzw_dW/m,
                           -w - dFzq_dq/m, -g*np.cos(theta)],
                          [q - dFyu_dU/m - dFytail_du/m, -dFyw_dW/m - dFytail_dw/m,
                           u - dFyq_dq/m - dFytail_dq/m, -g*np.sin(theta)],
                          [dMy_du/0.1, dMy_dw/0.1, dMy_dq/0.1, 0],
                          [0, 0, 1, 0]], float)
            return A

The main file
*************

The main is composed by the following parts:

1. Generate the kinematics object;
2. Build the **bird** object
3. Pass the bird object to the multiple-shooting solver
4. Solve the multiple-shooting problem

The main file is thus

.. code-block:: python

   
        import numpy as np
        from  ms_package.multiple_shooting import MultipleShooting
        from aero_package.bird_model import BirdModel
        from aero_package.kinematics_constructor import Shoulder, Elbow, Wrist, Joint
        import time
        from numpy.linalg import norm
        from numpy import inf
        from ms_package.lma_solver import Solver


        # generate bird kinematics by calling "kinematics_constructor" module
        bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                                 axis_y=Joint(-np.deg2rad(19),np.deg2rad(20),np.pi/2),
                                 axis_z=Joint(0,np.deg2rad(42),np.pi))
        bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                           axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

        bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                           axis_z=Joint(0.,0.,0.))

        mybird = BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist)

        # set initial guess for multiple-shooting scheme
        x0 = [18., 0.5, 0.1, 0.01]

        # generate multiple-shooting object
        ms_obj = MultipleShooting(x0, M = 2, period=0.25, t_steps=60, model = mybird)

        # call the LMA solver
        mysol = Solver(ms_obj = ms_obj).lma()

It is now possible to run the main file inside the ``multiflap`` directory::
  
   python3 main.py
