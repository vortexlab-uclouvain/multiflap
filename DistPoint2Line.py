import numpy as np

"""
Checked
"""

def DistPoint2Line(p, p0, v):
    s = np.sum(v*(p-p0))
    closest = p0 + (s*v)   # This might be .dot()
    d = np.sqrt(np.sum((p-closest)*(p-closest)))
    
    return s, closest, d