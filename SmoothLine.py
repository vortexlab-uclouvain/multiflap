import numpy as np  # Import NumPy

def SmoothLine(line, npt):
    
    smooth_radius = 0.05
    
    # Getting the length of the curve
    notok = 1
    l = np.zeros((3, npt))
    whilecicle = 0
    while notok == 1:
        whilecicle = whilecicle + 1
        n = np.size(line[0])
        lcurve = 0
        for i in range(n-1):
            lcurve = lcurve + np.sqrt(np.sum((line[:,i+1]-line[:,i])**2))
        
        ds = lcurve/(npt-1)
        
        l[:,0] = line [:,0]
        
        j = 0  #Until here is OK
        
        for i in range(1,npt):
            dist = np.sqrt(np.sum((line[:,j+1]-l[:,i-1])**2))
            while dist < ds:
                j = j+1
                if j < (np.size(line[0])-1):
                    dist = dist + np.sqrt(np.sum((line[:,j+1]-line[:,j])**2))
                else:
                    dist = ds + 1
            
            if j < (np.size(line[0])-1):
                fac = dist - ds
                den = np.sqrt(np.sum((line[:,j+1]-line[:,j])**2))
                fac = fac/den
                l[:,i] = (1-fac)*line[:,j+1]+fac*line[:,j]
            else:
                l[:,i] = line [:,-1]
        
        l2 = l
        
        notok = 0
        
        for i in range(1,npt-1):
            vec1 = l[:,i] - l[:,i-1]
            vec2 = l[:,i+1] - l[:,i]
            coss = np.sum(vec1*(vec2))/((np.linalg.norm(vec1))*(np.linalg.norm(vec2)))
            # =============================================================================
            # since truncation errors might lead to having coss slightly greater than 1
            # =============================================================================
            if coss > (1):
                angle = np.arccos(1.0)
            else:
                angle = np.arccos(coss)
            d = np.linalg.norm(vec1)
            smooth_angle = d/smooth_radius
            while angle > smooth_angle:
                notok = 1
                point = (l[:,i+1] + l[:,i-1])/2
                l2[:,i] = (l2[:,i] + point)/2
                vec1 = l2[:,i] - l[:, i-1]
                vec2 = l[:,i+1] - l2[:,i]
                coss = np.sum(vec1*(vec2))/((np.linalg.norm(vec1))*(np.linalg.norm(vec2)))
                angle = np.arccos(coss)
                d = np.linalg.norm(vec1)
                smooth_angle = d/smooth_radius
                
        line = l2
            
    return line
