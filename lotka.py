import numpy as np
from scipy.integrate import odeint
def lotka(x0, t):
    
    x, y= x0
    dxdt = 1.5*x - x*y
    dydt = x*y - 3*y

    vel_array = np.array([dxdt, dydt])
    return vel_array


time_array= np.linspace(0, 20, 1000)
intitial = np.array([1., 1.])
intitial_2 = np.array([5., 0.5])
ssp_soliution=odeint(lotka, intitial, time_array)
ssp_soliution_2 = odeint(lotka, intitial_2, time_array)
import matplotlib.pyplot as plt

fig1= plt.figure()
plt.plot(ssp_soliution[:,0],ssp_soliution[:,1], color='red')
plt.plot(ssp_soliution_2[:,0],ssp_soliution_2[:,1], color='blue')
ax1 = fig1.gca()
plt.xlabel('time (years)')
plt.ylabel('number of deers')
plt.grid()
plt.show()
