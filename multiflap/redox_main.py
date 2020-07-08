import numpy as np
from  odes.redox_oscillation import RedoxModel
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting_period import MultipleShooting
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ms_package.lma_solver_period import Solver

x = [0.5, 0.5, 0.6, 0.2]

time_array = np.linspace(0, 180, 90000)
mymodel = RedoxModel()

ms_obj =  MultipleShooting(x, M=2, period_guess= 23., t_steps=50000, model=mymodel)

mysol = Solver(ms_obj = ms_obj).lma(23., '/Users/gducci/UCL/PROJECT/Simulations/class_test')

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)


sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

plt.plot( sol_time, sol_array[:,0], label = "D1")
plt.plot( sol_time, sol_array[:,1], label = "D2")
plt.plot( sol_time, sol_array[:,2], label = "R")
plt.plot( sol_time, sol_array[:,3], label = "A")
plt.legend()
plt.show()
