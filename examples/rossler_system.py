import numpy as np
import multiflap as mf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
x = [.1, 5., 3.6]

time_array = np.linspace(0, 180, 90000)
mymodel = mf.Rossler(a=0.2, b=0.2, c=5.7)

ms_obj =  mf.MultipleShootingPeriod(x, M=30, period_guess= 5.,
                                    t_steps=5000, model=mymodel, integrator='odeint')

mysol = mf.SolverPeriod(ms_obj = ms_obj).lma()

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)

# plot the trajectory related to the initial condition
[_, initial_traj] = ms_obj.get_mappedpoint(x, 0, 100)
initial_points = ms_obj.get_initial_guess()

sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

initial_traj = initial_traj.x

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'b')
ax.scatter(initial_points[:, 0],
        initial_points[:, 1],
        initial_points[:, 2],color = 'red', alpha=0.5)
plt.legend()
plt.show()
fig2 = plt.figure(2)
ax = fig2.gca(projection='3d')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'b', label='limit cycle')
ax.plot(initial_traj[:, 0],
        initial_traj[:, 1],
        initial_traj[:, 2],color = 'red', alpha=0.5, label='initial point trajectory')
plt.legend()
plt.show()
