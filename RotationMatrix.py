import numpy as np  # Import NumPy

def RotationMatrix(theta, axis):
    
    if axis == "x":
        M = np.array([[1, 0, 0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]])
    elif axis == "y":
        M = np.array([[np.cos(theta), 0, np.sin(theta)],
                       [0, 1, 0],
                       [-np.sin(theta), 0, np.cos(theta)]])
    elif axis == "z":
         M = np.array([[np.cos(theta), -np.sin(theta), 0],
                        [np.sin(theta),  np.cos(theta),  0],
                        [0, 0, 1]])
    else:
        print("The input selected is not an axis")
    
    return M
