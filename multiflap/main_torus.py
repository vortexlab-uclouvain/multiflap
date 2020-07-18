import numpy as np
from  odes.torus import Torus
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting_period import MultipleShootingPeriod
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import SolverPeriod
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
# generate the ODEs object

mymodel = Torus(lam=2.)


## initial condition
x = [0.9, 0.9, 0.9]


# generate the multiple shooting object
#`ms_obj = MultipleShootingPeriod(x, M=3, period_guess=200.,
#`                                t_steps=10000, model=mymodel,
#`                                option_jacobian='numerical', integrator='odeint')
#`
#`
#`#guess = ms_obj.get_initial_guess()
#`#print("guess ok")
#`
#`# call the solver for the multiple-shooting algorithm
#`mysolution = SolverPeriod(ms_obj=ms_obj).lma()
#`#
#`jacobian = mysolution[4]
#`#
#`## Floquet multipliers
#`eigenvalues, eigenvectors = np.linalg.eig(jacobian)
#`## ODE limit cycle solution
#`sol_array = mysolution[3].space
#`sol_time = mysolution[3].time
#`period = sol_time[-1]
#sol_array = odeint(mymodel.dynamics, x, tarray)
# plot the phase portrait of the limit cycle
tarray=np.linspace(0,3000,100000)
sol_array= rk4(mymodel.dynamics, x, tarray)
sol_array = sol_array.x
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
#ax.scatter(initial_guess[:,0],
#           initial_guess[:,1],
#           initial_guess[:,2], color='red', label='initial guess')
ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'b', alpha=.3)
plt.legend()
plt.savefig('../img/isothermal_reaction.png')
plt.show()
