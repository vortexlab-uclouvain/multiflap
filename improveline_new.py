import numpy as np  # Import NumPy
from smooth_quarterline import smooth_quarterline
from distances_functions import distance_point2line, distance_line2line


def improveline_new(line, chord_leadingedge, chord_trailingedge, leadingedge, trailingedge):
    n = np.size(line[0])
    nle = np.size(leadingedge[0])
    nte = np.size(trailingedge[0])

    newline = np.zeros_like(line)
    newchordle = np.zeros_like(chord_leadingedge)
    newchordte = np.zeros_like(chord_trailingedge)
    newline[:,0] = line[:,0]
    newline[:,-1] = line[:,-1]
    newchordle[:,0] = chord_leadingedge[:,0]
    newchordle[:,-1] = chord_leadingedge[:,-1]
    newchordte[:,0] = chord_trailingedge[:,0]
    newchordte[:,-1] = chord_trailingedge[:,-1]

    # External Iteration

    for i in range(1,n-1):
        linevec = line[:,i+1] - line[:,i-1]
        linevec = linevec/np.linalg.norm(linevec)

        # orthogonalizing chordvec wrt linevec
        chordvec = chord_trailingedge[:,i] - chord_leadingedge[:,i]
        chordvec = chordvec - np.dot((np.sum(linevec*chordvec)),linevec)    # BUG HERE
        chordvec = chordvec/np.linalg.norm(chordvec)

        """
        getting points in le closest to chord line
        take the first two points of le as first references
        """
        jref1 = 0
        jref2 = 1

        [_ , _, dref1] = distance_point2line(leadingedge[:,jref1], line[:,i], chordvec) # If it does not work DP2L[][2] LAST ELEMENT OF OUTPUT
        [_, _, dref2] = distance_point2line(leadingedge[:,jref2], line[:,i], chordvec)


        """
        First IF Routine, checked (Below)
        """
        if (dref2<dref1):
            a = jref2
            jref2 = jref1
            jref1 = a
            a = dref2
            dref2 = dref1
            dref1 = a

        for j in range(2, nle):
            [_,_,dtest] = distance_point2line(leadingedge[:,j], line[:,i], chordvec)
            if (dtest < dref2):
                dref2 = dtest
                jref2 = j
            if (dref2 < dref1):
                a = jref2
                jref2 = jref1
                jref1 = a
                a = dref2
                dref2 = dref1
                dref1 = a

        levec = leadingedge[:,jref2] - leadingedge[:,jref1]
        levec = levec/np.linalg.norm(levec)
        [_, _, newchordle[:,i], _, _] = distance_line2line(leadingedge[:,jref1],levec,line[:,i],chordvec)

        # Repeat this for the trailing edge

        jref1 = 0
        jref2 = 1

        [_, _, dref1] = distance_point2line(trailingedge[:,jref1], line[:,i], chordvec) # If it does not work DP2L[][2] LAST ELEMENT OF OUTPUT
        [_, _, dref2] = distance_point2line(trailingedge[:,jref2], line[:,i], chordvec)

        if dref2 < dref1:
            a = jref2
            jref2 = jref1
            jref1 = a
            a = dref2
            dref2 = dref1
            dref1 = a

        for j in range(2, nte):

            [_,_,dtest] = distance_point2line(trailingedge[:,j], line[:,i], chordvec)
            if (dtest < dref2):
                dref2 = dtest
                jref2 = j
            if (dref2 < dref1):
                a = jref2
                jref2 = jref1
                jref1 = a
                a = dref2
                dref2 = dref1
                dref1 = a

        tevec = trailingedge[:,jref2] - trailingedge[:,jref1]
        tevec = tevec/np.linalg.norm(tevec)
        [_, _, newchordte[:,i], _, _] = distance_line2line(trailingedge[:,jref1],tevec,line[:,i],chordvec)
        newline[:,i] = .75*newchordle[:,i] + .25*newchordte[:,i]

    return(newline, newchordle, newchordte)
