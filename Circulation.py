import numpy as np
import polar as polar

def Circulation(y,z,yc,zc,chord,angle_of_attack,cl_alpha,velocity_module,gamma):
    
    gamma_new = np.zeros_like(gamma)
    vs = np.zeros([3,np.size(gamma)])
    relfac = 0.1

    gamma_filament = np.c_[gamma[:,0], gamma[:,1:] - gamma[0:,:-1], -gamma[0,-1]]
    
    dz = z[1:] - z[0:-1]
    dy = y[1:] - y[0:-1]
    ds = np.sqrt((dy**2)+(dz**2))
    dz = dz/ds
    dy = dy/ds
    down = np.array([dz,-dy,np.zeros_like(dy)])
    
    n = np.size(yc)
    
    for i in range(n):
        v = np.array([[0],[0],[0]])
        
        for j in range(n+1):
            
            dy_1 = yc[i] - y[j]
            dz_1 = zc[i] - z[j]
            
            v = v + gamma_filament[0,j]*np.array([[dz_1],[-dy_1],[0]])/(dy_1**2 + dz_1**2)
        
        v = (v/(4*np.pi))
        vdown = np.sum(v[:,0]*down[:,i])
        gamma_new[0,i] = 0.5*velocity_module[i]*chord[i]*cl_alpha*(angle_of_attack[i] - vdown/velocity_module[i])
        vs[:,i] = v[:,0]
    
    gamma_new = relfac*gamma_new + (1-relfac)*gamma
    
    return gamma_new, vs

def CirculationPolar(y,z,yc,zc,chord,angle_of_attack,cl_alpha,velocity_module,gamma):
    
    gamma_new = np.zeros_like(gamma)
    vs = np.zeros([3,np.size(gamma)])
    relfac = .1

    gamma_filament = np.c_[gamma[:,0], gamma[:,1:] - gamma[0:,:-1], -gamma[0,-1]]
    
    dz = z[1:] - z[0:-1]
    dy = y[1:] - y[0:-1]
    ds = np.sqrt((dy**2)+(dz**2))
    dz = dz/ds
    dy = dy/ds
    down = np.array([dz,-dy,np.zeros_like(dy)])
    
    n = np.size(yc)
    
    for i in range(n):
        v = np.array([[0],[0],[0]])
        
        for j in range(n+1):
            
            dy_1 = yc[i] - y[j]
            dz_1 = zc[i] - z[j]
            
            v = v + gamma_filament[0,j]*np.array([[dz_1],[-dy_1],[0]])/(dy_1**2 + dz_1**2)
        
        v = (v/(4*np.pi))
        vdown = np.sum(v[:,0]*down[:,i])
        gamma_new[0,i] = 0.5*velocity_module[i]*chord[i]*polar.lift_coefficient(np.rad2deg(angle_of_attack[i] - vdown/velocity_module[i]))
        vs[:,i] = v[:,0]
    
    gamma_new = relfac*gamma_new + (1-relfac)*gamma
    
    return gamma_new, vs

def InducedVel_Filament(X1,Xleg,legdir, vstrength):
### returns the induced velocity on X1 of a seminfinite circulation 
###leg starting in Xleg directed along legdir and with strenght vstrength (BIOT-Savart law integrated from Xleg to \infty)
### for using in lifting line legdir is direction of Vinf, vstrength is d\lambda
    tol=1e-12
    legdir=legdir/np.linalg.norm(legdir)
    
    R1leg=X1-Xleg
    Rmod=np.linalg.norm(R1leg)
    
    rdir=R1leg/Rmod
    crossp=np.cross(legdir,rdir)/np.linalg.norm(np.cross(legdir,rdir))

    d=np.linalg.norm(R1leg-np.dot(R1leg,legdir)*legdir)
    if np.abs(d)<tol: vind=np.zeros(crossp.shape)
    else:
        V1coeff=1/4/np.pi/d*(1+np.dot(rdir,legdir))*crossp
        vind=vstrength*V1coeff
    return(vind)

def InducedVel_Segment(X,Xlegs,Xlege,vstrength,sigma):
### calculates the exact regularized version given by Wincklemans (LOA Rosenhead-Moore)
    ### X is where you calculate the induced velocity. Xlegs/e are starting and ending point of the segment, vstrength is the circulation.
    ### sigma is a regularization parameter which by default is zero
    ## returns induced velocity 3D vind
    ## sigma == 0
    
    legdir=(Xlege-Xlegs)/np.linalg.norm(Xlege-Xlegs)
    rs=X-Xlegs
    h=np.linalg.norm(rs-np.dot(rs,legdir)*legdir)
    z2=np.dot((X-Xlege),legdir)
    z1=np.dot((X-Xlegs),legdir)
    direc=np.cross(legdir,rs); dist=np.linalg.norm(direc)
    vind=np.zeros(3)
    if dist>1e-12:
        vdir=direc/dist
        Start=z1/(z1**2+h**2+sigma**2)**.5
        End=z2/(z2**2+h**2+sigma**2)**.5
        vcoeff=-vstrength/4/np.pi*h/(h**2+sigma**2)*(End-Start)
        vind=vcoeff*vdir
    return(vind)