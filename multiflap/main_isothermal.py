import numpy as np
from  odes.isothermal_reaction import IsothermalReaction
from ms_package.multiple_shooting_period import MultipleShootingPeriod
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import SolverPeriod
from mpl_toolkits.mplot3d import Axes3D

# generate the ODEs object

mymodel = IsothermalReaction(lam=11.)


# initial condition
x = [40., 20., 20.]

# generate the multiple shooting object
ms_obj = MultipleShootingPeriod(x, M=20, period_guess=.5,
                                t_steps=200, model=mymodel)

# just to plot the initial guess distribution. No need to call this
initial_guess = ms_obj.get_initial_guess()

# call the solver for the multiple-shooting algorithm
mysolution = SolverPeriod(ms_obj=ms_obj).lma()

jacobian = mysolution[4]

# Floquet multipliers
eigenvalues, eigenvectors = np.linalg.eig(jacobian)

# ODE limit cycle solution
sol_array = mysolution[3].space
sol_time = mysolution[3].time
period = sol_time[-1]

# plot the phase portrait of the limit cycle
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
ax.scatter(initial_guess[:,0],
           initial_guess[:,1],
           initial_guess[:,2], color='red', label='initial guess')
ax.plot(sol_array[:, 0],
        sol_array[:, 1],
        sol_array[:, 2],color = 'b')
plt.legend()
plt.show()
