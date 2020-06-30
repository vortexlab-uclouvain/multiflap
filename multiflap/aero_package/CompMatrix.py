import numpy as np  # Import NumPy

def CompMatrix(k, axis):
    
    if axis == "x":
        M = np.array([[k, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
    elif axis == "y":
        M = np.array([[1, 0, 0],
                      [0, k, 0],
                      [0, 0, 1]])
    elif axis == "z":
          M = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, k]])
    else:
        print("The input selected is not an axis")
    
    return M
