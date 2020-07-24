'''
file distances_functions.py
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
import numpy as np

def distance_point2line(p, p0, v):
    s = np.sum(v*(p-p0))
    closest = p0 + (s*v)   # This might be .dot()
    d = np.sqrt(np.sum((p-closest)*(p-closest)))

    return s, closest, d



def distance_line2line(p01, v1, p02, v2):

    v1v2 = np.sum(v1*(v2))
    p01v1 = np.sum(p01*(v1))
    p02v1 = np.sum(p02*(v1))
    p01v2 = np.sum(p01*(v2))
    p02v2 = np.sum(p02*(v2))

    s1 = (p02v1 - p01v1 + p01v2*v1v2 - p02v2*v1v2)/(1 - v1v2*v1v2)
    s2 = p01v2 + s1*v1v2 - p02v2

    p1 = p01 + s1*v1
    p2 = p02 + s2*v2
    d = np.sqrt(np.sum((p1-p2)*((p1-p2))))

    return(s1, s2, p1, p2, d)
