'''
file flapping_forces.py
@author Victor Colognesi, Gianmarco Ducci
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
from .MergeLines import MergeLines
from .circulation import circulation, inucedvel_filament, inducedvel_segment
from .TailModels import delta_tail, tail_geometry
from .settings import SimulationsSettings
from .bird import Shoulder, Elbow, Wrist, Joint
settings = SimulationsSettings()
def get_flappingforces(x0, v_kin, line, chordir, updir, chord, tail_opening):

    u, w, q, theta = x0
    cl_alpha = 2*np.pi
    rho = 1.225  # Air density
    CD_0 = 0.02
    F_pro = 0
    U = np.array([0, w, u])
    ang_velocity = np.array([q, 0, 0])

    wingframe_position = settings.wingframe_position
    tail_frame_position = settings.wingframe_position_tail
    # -------------------------------------------------------------------------
    # Initialization of the vector force. It must have the same dimension of the time vector
    # -------------------------------------------------------------------------
    Fx = 0
    Fy = 0
    Fz = 0
    My = 0
    M_wing = 0
    M_drag = 0
    M_lift = 0
    M_tail = 0
    F_tail_tot = 0
    P_ind = 0

    line_c = line[:, 1:-1:2] # Center points: at the center of each
    # interval --> every other point,
    # starting from the second until the one
    # before the last one
    updir = updir[:, 1:-1:2]  # updir, chordir and chord at the center points
    chordir = chordir[:, 1:-1:2]
    chord = chord[1:-1:2]
    line = line[:, ::2]  # "mid" points: points at the junction between two segments --> every other point from the first one to the last one

    # -------------------------------------------------------------------------
    #  u_inf = free stream velocity, equal over each profile (points)
    # -------------------------------------------------------------------------

    u_inf = np.array([U for j in range(np.size(line_c[0]))]).T
    # tangential_wing = (line_c - line_c2)/dt
    line_velocity = np.zeros_like(line_c)
    for j in range(np.size(line_c[0])):
        line_velocity[:, j] = np.cross(ang_velocity, wingframe_position+line_c[:, j])

    line_velocity = line_velocity + v_kin

#    velocity = u_inf - (line_c - line_c2)/dt

    velocity = u_inf - line_velocity

    velocity_profile = np.zeros_like(velocity)
    line_direction = np.zeros_like(updir)
    norm_velocity = np.zeros_like(chord)
    direction_velocity = np.zeros_like(velocity)
    angle_of_attack = np.zeros_like(chord)

    for j in range(np.size(line_c[0])):
        line_direction[:, j] = np.cross(updir[:, j], chordir[:, j])
        velocity_profile[:, j] = velocity[:, j] - line_direction[:, j]*np.sum(line_direction[:, j]*velocity[:, j])
        norm_velocity[j] = np.sqrt(np.sum(velocity_profile[:, j]*velocity_profile[:, j]))
        direction_velocity[:, j] = velocity_profile[:, j]/norm_velocity[j]
        angle_of_attack[j] = np.arctan2(np.sum(direction_velocity[:, j]*updir[:, j]), np.sum(direction_velocity[:, j]*chordir[:, j]))
    tol = 0.005
    max_iterations = 50
    ds = np.zeros(len(chord))
    s = np.zeros((np.size(line[0])))
#    print(np.mean(np.rad2deg(angle_of_attack)))
    # -------------------------------------------------------------------------
    # Evaluation of ds length by distance point to point
    # ds is the distance between two consecutives discretised points
    # that belong to the lifting line
    # -------------------------------------------------------------------------

    for j in range(np.size(line[0])-1):
        ds[j] = np.sqrt(np.sum((line[:, j+1] - line[:, j])**2))

    for j in range(1, np.size(line[0])):
        s[j] = s[j-1] + ds[j-1]

    gamma = np.zeros((1, np.size(ds)))
    error = tol + 1
    error_p = np.inf
    count = 0

    while (error > tol) and (count < max_iterations):

        count = count + 1

        [gamma_new, v_downwash] = circulation(line[0,:], line[1,:], line_c[0,:], line_c[1,:],chord, angle_of_attack, cl_alpha, norm_velocity, gamma)

        error = np.mean(np.abs(gamma - gamma_new))/np.mean(np.abs(gamma_new))

        if error > error_p and count > 0:
            print("Error increasing after iteration ",count)
            break

        error_p = error

        gamma = gamma_new

    gamma_filament = np.c_[gamma[:,0], gamma[:,1:] - gamma[0:,:-1], -gamma[0,-1]]
    gamma = np.reshape(gamma, np.size(gamma))
    gamma_filament = np.reshape(gamma_filament, np.size(gamma_filament))
    local_velocity = velocity + v_downwash
    lever_arm = np.zeros_like(line_c)
    BS_radius = np.zeros_like(line_c)
# =============================================================================
#     Induced velocity on the tail
# =============================================================================
    [tail_span, AR_tail, NP_tail] = tail_geometry(tail_opening)
    V_ind = 0
    dl_induced = line[:,1:] - line[:,:-1]
    for i in range (np.size(gamma)):
        lever_arm[:, i] = line_c[:, i] + wingframe_position
        BS_radius[:,i] = (NP_tail - lever_arm[:, i])
        cross_prod = (np.cross(dl_induced[:,i],BS_radius[:,i]))
        V_ind = V_ind + ((gamma[i]/(4*np.pi))*cross_prod)/(np.linalg.norm(BS_radius[:,i])**3)
#        V_ind = V_ind + circulation.InducedVel_Segment(NP_tail,line[:,i]+wingframe_position, line[:,i+1]+wingframe_position, gamma[i], 1e-13)

    lever_arm1 = np.zeros_like(line)
    BS_radius1 = np.zeros_like(line)
    for i in range (len(line)):
        lever_arm1[:, i] = line[:, i] + wingframe_position
        BS_radius1[:,i] = (NP_tail - lever_arm1[:, i])
        V_ind = V_ind + inucedvel_filament(BS_radius1[:,i],line[:,i], U, gamma_filament[i])


    velocity_tail_q = np.cross(ang_velocity,tail_frame_position)
    U_tail = U + V_ind + velocity_tail_q
#    U_tail = U + V_ind
    M = np.zeros_like(gamma)
#    alpha_tail = np.arctan2(U_tail[1],U_tail[2])
    for j in range(np.size(gamma)):
        F = rho*gamma[j]*np.cross(local_velocity[:,j] , line_direction[:,j])
        lever_arm[:, j] = line_c[:, j] + wingframe_position
        F_moment = np.array([F[0]*ds[j], F[1]*ds[j], F[2]*ds[j]])
        Drag_moment = np.array([0., 0., F[2]*ds[j]])
        Lift_moment = np.array([0., F[1]*ds[j], 0.])

        Fx = Fx + F[0]*ds[j]  # Component along wing span
        Fy = Fy + F[1]*ds[j]  # Component of Lift (perpendicular to the free stream)
        Fz = Fz + F[2]*ds[j]  # Component of drag (along free stream)
        M_drag_j = np.cross(lever_arm[:, j], Drag_moment)
        M_lift_j = np.cross(lever_arm[:, j], Lift_moment)
        M_wing_j = np.cross(lever_arm[:, j], F_moment)
        F_pro = F_pro + (0.5*rho*(CD_0)*chord[j]*local_velocity[:,j]**(2))*ds[j]
        M_wing = M_wing + M_wing_j[0]
        M_drag = M_drag + M_drag_j[0]
        M_lift = M_lift + M_lift_j[0]
        My = My + M[0]
        P_ind = P_ind + np.dot(F*ds[j],-line_velocity[:,j])
        #P_ind = P_ind + np.dot(F*ds[j],U)

    F_tail_tot = delta_tail(U_tail, tail_span)
    M_tail = np.cross(NP_tail, F_tail_tot)
    My = M_wing + M_tail[0]
    return Fx, Fy, Fz+F_pro[2], My, F_tail_tot[1], F_tail_tot[2], M_wing, M_tail[0], M_drag, M_lift, P_ind, F_pro[2], angle_of_attack
