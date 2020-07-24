'''
file MergeLines.py
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

def MergeLines(line_r,updir_r,chordir_r,chord_r,line_l,updir_l,chordir_l,chord_l,nmid,x_ep):
    
    line_r[0,:] = line_r[0,:] + x_ep 
    line_l[0,:] = line_l[0,:] - x_ep 
    
    xmid = np.linspace(line_l[0,-1],line_r[0,0],nmid)
    length_xmid = np.size(xmid)
    line_mid = np.array([xmid, line_r[1,0]*np.ones(np.size(xmid)), line_r[2,0]*np.ones(np.size(xmid))])
    chordir_mid = np.array([chordir_r[:,0] for i in range(np.size(xmid))]).T
    updir_mid = np.array([updir_r[:,0] for i in range(np.size(xmid))]).T
    chord_mid = np.zeros(length_xmid)
    chord_mid[:] = chord_r[0]*np.ones(length_xmid)
    
    line = np.c_[line_l,line_mid[:,1:-1],line_r]
    chordir = np.c_[chordir_l, chordir_mid[:,1:-1], chordir_r]
    updir = np.c_[updir_l, updir_mid[:,1:-1], updir_r]
    chord = np.concatenate([chord_l, chord_mid[1:-1]/1e4, chord_r])
    return line, chordir, updir, chord
