import numpy as np
import multiflap as mf
import matplotlib.pyplot as plt
from scipy.integrate import odeint
x = [1., 1.]
my_model = mf.ForcedVanDerPol(lam=1.8)

ms_obj = mf.MultipleShooting(x, M=2, period=10.4719755, t_steps=4000,
                             model=my_model, option_jacobian='numerical')

mysolution = mf.Solver(tolerance=1e-6, ms_obj = ms_obj).lma()

sol_array = mysolution[3].space
sol_time = mysolution[3].time
period = sol_time[-1]

x0 = sol_array[0,:]
jac = mysolution[4]
eigenvalues, eigenvectors = np.linalg.eig(jac)

print(eigenvalues)
tarray = np.linspace(0, 10, 40000)
solution_odeint = odeint(my_model.dynamics, x0, tarray)

plt.plot( sol_array[:,0], sol_array[:,1])
plt.plot(solution_odeint[:,0], solution_odeint[:,1], color='red')
plt.scatter(sol_array[0,0], sol_array[0,1])
plt.show()
