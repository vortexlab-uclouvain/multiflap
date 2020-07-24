'''
file CompMatrix.py
@author Victor Colognesi
@copyright Copyright © UCLouvain 2020

multiflap is a Python tool for finding periodic orbits and assess their stability via the Floquet multipliers.

Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>

List of the contributors to the development of multiflap, Description and complete License: see LICENSE and NOTICE files.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
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
