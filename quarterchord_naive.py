import numpy as np  # Import NumPy

def quarterchord_naive(leadingedge, trailingedge):
    line = .75*leadingedge + .25*trailingedge;
    chord_leadingedge = leadingedge;
    chord_trailingedge = trailingedge;
    return line, chord_leadingedge, chord_trailingedge
