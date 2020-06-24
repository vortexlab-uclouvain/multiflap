import numpy as np

"""
Checked
"""

def distance_point2line(p, p0, v):
    s = np.sum(v*(p-p0))
    closest = p0 + (s*v)   # This might be .dot()
    d = np.sqrt(np.sum((p-closest)*(p-closest)))

    return s, closest, d
