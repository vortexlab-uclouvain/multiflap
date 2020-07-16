import numpy as np
import multiflap as mf
import matplotlib.pyplot as plt

x = [10., 10., 3.6]

time_array = np.linspace(0, 180, 90000)
mymodel = mf.Rossler(a=0.2, b=0.2, c=5.7)

ms_obj =  mf.MultipleShootingPeriod(x, M=2, period_guess= 5., t_steps=50000, model=mymodel)

mysol = mf.SolverPeriod(ms_obj = ms_obj).lma(5.)

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)


sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

plt.plot( sol_time, sol_array[:,0], label = "D1")
plt.plot( sol_time, sol_array[:,1], label = "D2")
plt.plot( sol_time, sol_array[:,2], label = "R")
plt.legend()
plt.show()
