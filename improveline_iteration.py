import numpy as np  # Import NumPy
from smooth_quarterline import smooth_quarterline
from improveline_new import improveline_new

def improveline_iteration(line, chord_leadingedge, chord_trailingedge, leadingedge, trailingedge, it):

    relfactor = .25

    [line, chord_leadingedge, chord_trailingedge] = smooth_quarterline(line,chord_leadingedge,chord_trailingedge)
    [newline, newchordle, newchordte] = improveline_new(line, chord_leadingedge, chord_trailingedge,leadingedge,trailingedge)

    chgt = 0

    for i in range(np.size(line[0]) - 1):

        local_max = np.maximum(0, i-1)
        local_min = np.minimum(np.size(line[0]), i + 1)
        norm1 = np.linalg.norm(newline[:,i] - line[:,i])
        norm2 = np.linalg.norm(line[:,local_max]-line[:,local_min])
        chgt = chgt + relfactor*norm1/norm2

    chgt = chgt/np.size(line[0])

    line = relfactor*newline + (1-relfactor)*line
    chord_leadingedge = relfactor*newchordle + (1-relfactor)*chord_leadingedge
    chord_trailingedge = relfactor*newchordte + (1-relfactor)*chord_trailingedge

    return line, chord_leadingedge, chord_trailingedge, chgt
