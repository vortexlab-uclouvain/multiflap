import numpy as np
from  ms_package.redox_oscillation import RedoxModel
from ms_package.rk_integrator import rk4
from ms_package.multiple_shooting import MultipleShooting
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ms_package.lma_solver import Solver

x = [0.1, 0.5, 0.6, 0.2]

time_array = np.linspace(0, 180, 90000)
mymodel = RedoxModel()

ms_obj =  MultipleShooting(x, M=4, period = 24.952, t_steps=50000, model = mymodel)

mysol = Solver(ms_obj = ms_obj).lma('/Users/gducci/UCL/PROJECT/Simulations/class_test')

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)


sol_array = mysol[3].space
sol_time = mysol[3].time
plt.plot( sol_time, sol_array[:,0])
plt.plot( sol_time, sol_array[:,1])
plt.plot( sol_time, sol_array[:,2])
plt.plot( sol_time, sol_array[:,3])

plt.show()