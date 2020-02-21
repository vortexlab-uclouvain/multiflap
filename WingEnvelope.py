import numpy as np  # Import NumPy
from RotationMatrix import RotationMatrix
from CompMatrix import CompMatrix
from Plumes import plumes


def WingEnvelope(shoulder_x, shoulder_z, shoulder_y, elbow_y, elbow_x, wrist_y, wrist_z):
    
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
       # Ellipse thing (just double check to delete unnecessary things)
    number_ellipses = 2
    s = np.linspace(0,2*np.pi,number_ellipses + 1) # No idea but works
    s = s[0:number_ellipses]
    skz = .01
    sky = .01 
    z = skz*np.cos(s)
    y = sky*np.sin(s)
    points = np.zeros((3 , nslice*number_ellipses));
    
    # Weight thing, again ask to have a more clear comment
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
                faccd1 = 0.5*(1-np.sin(np.pi*(xslice[j]-l1)/lrefcd));
         
            faccd2 = 1-faccd1;
            fac2a = 1;
            fac2b = 0;
            w1[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd1
            w2a[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2a
            w2b[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2b
            w3[j*number_ellipses:(j+1)*number_ellipses] = (1-facb)
            w4[j*number_ellipses:(j+1)*number_ellipses] = 0
            
        elif xslice[j] < l1 + l2:
            for k in range(number_ellipses):
                points[:,(j)*number_ellipses+k] = [xslice[j], y[k], z[k]]

            facb = (z/skz + 1)/2;
            lrefcd = min(l1,l2);
            if xslice[j] < l1+lrefcd/2:
                faccd2 = 0.5*(1+np.sin(np.pi*(xslice[j]-l1)/lrefcd));
            else:
                faccd2 = 1;
        
            faccd1 = 1-faccd2;
            fac2b = (xslice[j]-l1)/l2;
            fac2a = 1-fac2b;
            lrefpg = min(l2,l4);
            
            if xslice[j] < l1+l2-lrefpg/2:
                facpg2 = 1;
            else:
                facpg2 = 0.5*(1-np.sin(np.pi*(xslice[j]-(l1+l2))/lrefpg));
        
            facpg4 = 1-facpg2;
        
            w1[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd1*facpg2
            w2a[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2a*facpg2
            w2b[j*number_ellipses:(j+1)*number_ellipses] = facb*faccd2*fac2b*facpg2
            w3[j*number_ellipses:(j+1)*number_ellipses] = (1-facb)*facpg2
            w4[j*number_ellipses:(j+1)*number_ellipses] = facpg4
        else:
            for k in range(number_ellipses):
                points[:,(j)*number_ellipses+k] = [xslice[j], y[k], z[k]]
        
            facb = (z/skz + 1)/2;
            fac2a = 0;
            fac2b = 1;
            lrefpg = min(l2,l4);
            
            if xslice[j] < l1+l2+lrefpg/2:
                facpg4 = 0.5*(1+np.sin(np.pi*(xslice[j]-(l1+l2))/lrefpg))
            else:
                facpg4 = 1;
            
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
            prerot = rot4;
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
            prerot = rot1;
    
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


if __name__ == "__main__":
   [LE, TE] = WingEnvelope(1, 1, 1, 0, 0, 0, 5)
