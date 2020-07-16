import numpy as np
import multiflap as mf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = [10., 10., 3.6]

mymodel = mf.Lorentz(a=10, b=28, c=8/3)

ms_obj =  mf.MultipleShootingPeriod(x, M=2, period_guess= 10., t_steps=50000, model=mymodel)

odes = mymodel.dynamics

t_array = np.linspace(0., 50, 50000)
sol = mf.rk4(odes, x, t_array)

mysol = mf.SolverPeriod(ms_obj = ms_obj).lma(5.)

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)

sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

sol_integration = sol.x

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
ax.plot(sol_integration[:, 0],
        sol_integration[:, 1],
        sol_integration[:, 2], alpha = 0.3, color='red')

ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'b')
plt.show()
