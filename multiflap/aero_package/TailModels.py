'''
file TailModels.py
@author Gianmarco Ducci
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
from .settings import SimulationsSettings

settings = SimulationsSettings()
rho = 1.225
kin_viscosity = 1.5e-5
Reynolds = float(15*settings.tail_length/kin_viscosity)

def delta_tail(velocity_vector, tail_span):
    alpha = np.arctan(velocity_vector[1]/velocity_vector[2])
    velocity_module = np.linalg.norm(velocity_vector)
    Force = (np.pi/4)*rho*(velocity_module**2)*alpha*(tail_span**2)
    Lift = Force*np.cos(alpha)
    Cf = 1.328/np.sqrt(Reynolds)
    Df = 0.5*rho*(velocity_module**2)*(tail_span*settings.tail_length/2)*Cf
    Drag_ind = 0.5*Lift*alpha
    Drag = Drag_ind + Df + Force*np.sin(alpha)
    Lat_force = 0.
    return Lat_force, Lift, Drag

def tail_geometry(length, **opening):
    tail_opening=opening.get('tail_opening', settings.tail_opening)
    tail_span = 2*length*np.tan(tail_opening/2)
    AR_tail = 2*(tail_span/length)
    NP_tail = settings.wingframe_position_tail + np.array([0., 0., (2/3)*length])
    return tail_span, AR_tail, NP_tail


if __name__ == "__main__":
    [b, _, _] = tail_geometry(.25, tail_opening=np.deg2rad(40))
    print(b)
