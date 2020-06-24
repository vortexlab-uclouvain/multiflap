import numpy as np
"""
Cheched,works properly 21/02/2019
"""


def smooth_quarterline(oldline,oldchordle,oldchordte):

    npt = np.size((oldline[0]))
    lcurve = 0;
    for i in range(npt-1):
        lcurve = lcurve + np.sqrt(np.sum((oldline[:,i+1]-oldline[:,i])**2))

    ds = lcurve/(npt-1)

    line = np.zeros_like(oldline)
    chord_leadingedge = np.zeros_like((oldline))
    chord_trailingedge = np.zeros_like((oldline))
    line[:,0] = 1*oldline[:,0]
    chord_leadingedge[:,0] = oldchordle[:,0]
    chord_trailingedge[:,0] = oldchordte[:,0]
    j = 0

    for i in range(1,npt):
        dist = np.sqrt(np.sum((oldline[:,j+1]-line[:,i-1])**2))
        while dist < ds:
            j = j+1
            if j < (np.size(oldline[0])-1):
                dist = dist + np.sqrt(np.sum((oldline[:,j+1]-oldline[:,j])**2))
            else:
                dist = ds + 1

        if j < (np.size(line[0])-1):
            fac = dist - ds
            den = np.sqrt(np.sum((oldline[:,j+1]-oldline[:,j])**2))
            fac = fac/den
            line[:,i] = (1-fac)*oldline[:,j+1]+fac*oldline[:,j]
            chord_leadingedge[:,i] = (1-fac)*oldchordle[:,j+1]+fac*oldchordle[:,j]
            chord_trailingedge[:,i] = (1-fac)*oldchordte[:,j+1]+fac*oldchordte[:,j]
        else:
            line[:,i] = 1*oldline [:,-1]
            chord_leadingedge[:,i] = oldchordle[:,-1]
            chord_trailingedge[:,i] = oldchordte[:,-1]

    return line, chord_leadingedge, chord_trailingedge

if __name__ == "__main__":
   smooth_quarterline(1,1,1)
