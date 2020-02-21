import numpy as np

def DistLine2Line(p01, v1, p02, v2):
    
    v1v2 = np.sum(v1*(v2)) 
    p01v1 = np.sum(p01*(v1)) 
    p02v1 = np.sum(p02*(v1)) 
    p01v2 = np.sum(p01*(v2)) 
    p02v2 = np.sum(p02*(v2)) 
    
    s1 = (p02v1 - p01v1 + p01v2*v1v2 - p02v2*v1v2)/(1 - v1v2*v1v2) 
    s2 = p01v2 + s1*v1v2 - p02v2 
    
    p1 = p01 + s1*v1 
    p2 = p02 + s2*v2 
    d = np.sqrt(np.sum((p1-p2)*((p1-p2))))
    
    
    
    
    return(s1, s2, p1, p2, d)
