"""
Phase space animation scritp

Script below is a template for creating an animation in the state-space domain
of a trjectory obtained by numerical integration. Rossler's system has been
encoded highlighting the unstable behaviour.

Copyright (C) <2019> <Université catholique de Louvain (UCL), Belgique>

The contributor to the development of this code is Gianmarco Ducci.

Website: https://sites.uclouvain.be/RevealFlight/gducci
Contact: gianmarco.ducci@uclouvain.be

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*************************************************************************
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from scipy.integrate import odeint

def Lorentz(initial_value, t):
    a = 10
    b = 28
    c = 8/3

    x, y, z = initial_value
    dxdt = a*(y-x)
    dydt = x*(b-z)-y
    dzdt = x*y - c*z
    eqns = np.array([dxdt, dydt, dzdt], float)

    return eqns


fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1,projection="3d")
ax1.view_init( 10, -65)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 4)


initial_value_lorentz = np.array([-2.39918211, -4.38992764, 10.22727544])
initial_value_lorentz_perturbed = initial_value_lorentz[:] + initial_value_lorentz[:]*10e-3
timestep = 4000
tArray = np.linspace(0, 30, timestep) # Time integration discretisation
test = odeint(Lorentz, initial_value_lorentz,tArray)  # Numerical Integration of the limit cycle unperturbed
perturbed_orbit = odeint(Lorentz, initial_value_lorentz_perturbed, tArray) #Numerical Integration of the limit cycle unperturbed

x = test[:,0]
y = test[:,1]
z = test[:,2]

x_perturbed = perturbed_orbit[:,0]
y_perturbed = perturbed_orbit[:,1]
z_perturbed = perturbed_orbit[:,2]


ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')

ax2.set_facecolor((1.0, 1.0, 0.95, 1.0))
ax2.set_xlabel('t')
ax2.set_ylabel('y')
ax2.set_xlim(0, tArray[-1])

ax3.set_facecolor((1.0, 1.0, 0.95, 1.0))
ax3.set_xlabel('t')
ax3.set_ylabel('z')
ax3.set_xlim(0, tArray[-1])

lines = []
for i in range(len(tArray)):
    head = i
    head_slice = (tArray > tArray[i] - .1) & (tArray< tArray[i])
    line1,  = ax1.plot(x[:i], y[:i], z[:i],
                       color='blue', alpha=.9)
    line1_inst,  = ax1.plot(x_perturbed[:i+1], y_perturbed[:i+1], z_perturbed[:i+1],
                       color='red', alpha=.5)
    line1_slice, = ax1.plot(x[head_slice], y[head_slice], z[head_slice],
                       color='blue', linewidth=2, alpha=.9)
    line1_head, = ax1.plot([x[head]], [y[head]], [z[head]],
                       color='blue', marker='o', markeredgecolor='blue')
    line1_head_inst, = ax1.plot([x_perturbed[head]], [y_perturbed[head]], [z_perturbed[head]],
                       color='red', marker='o', markeredgecolor='red')
    line2,  = ax2.plot(tArray[:i], y[:i],
                       color='blue', alpha=.9)
    line2_inst,  = ax2.plot(tArray[:i], y_perturbed[:i],
                       color='red', alpha=.5)

    line2_head, = ax2.plot(tArray[i], y[i],
                       color='blue', marker='o', markeredgecolor='blue')
    line2_head_inst, = ax2.plot(tArray[i], y_perturbed[i],
                       color='red', marker='o', markeredgecolor='red')
    line3,  = ax3.plot(tArray[:i], z[:i],
                       color='blue',  alpha=.9)
    line3_inst, = ax3.plot(tArray[:i], z_perturbed[:i],
                       color='red', alpha=.5)
    line3_head, = ax3.plot(tArray[i], z[i],
                       color='blue', marker='o', markeredgecolor='blue')
    line3_inst_head, = ax3.plot(tArray[i], z_perturbed[i],
                   color='red', marker='o', markeredgecolor='red')

    lines.append([line1,line1_inst,line1_head_inst, line1_slice,line1_head,line2,
                  line2_inst,line2_head,line2_head_inst, line3,line3_inst,line3_head,line3_inst_head])

plt.tight_layout()
ani = animation.ArtistAnimation(fig, lines, interval=50, blit=True)
ani.save('../../anim/phase_animation_3D.mp4', writer='ffmpeg',fps=4000/50)
plt.rcParams['animation.html'] = 'html5'
ani
