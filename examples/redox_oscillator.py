import numpy as np
import multiflap as mf
import matplotlib.pyplot as plt

x = [0.5, 0.5, 0.6, 0.2]

time_array = np.linspace(0, 180, 90000)
mymodel = mf.RedoxModel()

ms_obj =  mf.MultipleShootingPeriod(x, M=2, period_guess= 23.,
                                    t_steps=50000, model=mymodel)

mysol = mf.SolverPeriod(ms_obj = ms_obj).lma(23.)

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
