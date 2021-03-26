import numpy as np
import re
import glob
import os
import matplotlib.pyplot as plt
import matplotlib as mpl


def list_simulations(pathname):
    """
    pathname = /mypath/element_name
    example: pathname = './tail_*' to list all the elements starting for "tail_"
    """
    simlist = glob.glob(pathname)
    return simlist

sim = list_simulations('./tail_*_sweep*')

# Dimension an array of 0 and 1 flagging the cycles that are found
# [0/1, tail_opening, sweep_offset]
conditional_array = np.zeros((len(sim), 3))
list_velocities = []
list_multipliers = []
for i in range(len(sim)):
    simname = sim[i]
    # check wheter the limit cycle was found (this is rough!)
    # should put also a filter on the vertical velocity
    if os.path.isfile(simname+ '/solution_LC.npy'):
        x = re.findall(r'\d*\.?\d+', simname)
        x = list(map(float, x))
        u_fix = np.load(simname+ '/u_fix.npy')
        w_fix = np.load(simname+ '/w_fix.npy')
        solution = np.load(simname+'/solution_LC.npy' )
        eigenvalues = np.load(simname+'/floquet_multipliers.npy')
        theta = solution[:,-1]
        if np.abs(np.mean(w_fix)) < 0.05 and np.mean(u_fix) > 13.5 and np.abs(np.mean(theta)) < 0.15:
            conditional_array[i, 0] = 1
            conditional_array[i, 1:] = x
            solution_LC = np.load(simname+ '/solution_LC.npy')
            list_velocities.append(np.mean(u_fix))
            list_multipliers.append(eigenvalues)           
            print("u_f =", np.mean(u_fix))
        else:
            conditional_array[i, 1:] = x
    else:
        x = re.findall(r'\d*\.?\d+', simname)
        x = list(map(float, x))
        conditional_array[i, 1:] = x
        print(simname+": LC not detected")
# indices for 'bad' values
velocities = np.asarray(list_velocities)
colors = velocities**2
indices = conditional_array[:,0] == 0.
print(len(velocities))
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
x = conditional_array[~indices,1] # tail array
y = conditional_array[~indices,2] # sweep array
z = velocities
#ax1.tricontourf(x, y, z, levels=40,linewidths=0.5, colors='k')
#cntr2 = ax1.tricontourf(x, y, z, levels=40, cmap="RdBu_r")
cm = plt.cm.get_cmap('RdYlBu')
#fig1.colorbar(cntr2, ax=ax1)

sc = ax1.scatter(conditional_array[~indices,1], conditional_array[~indices,2], marker='s',  s=4,c=velocities, cmap = cm)
#ax1.plot(x, max(y))
plt.colorbar(sc)
#ax1.scatter(conditional_array[indices,1], conditional_array[indices,2], c='red')
ax1.set_ylabel('Sweep offset [$\deg$]')
ax1.set_xlabel('Tail opening [$\deg$]')
ax1.grid(False)
plt.show()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
x = conditional_array[~indices,1]
y = conditional_array[~indices,2]
z = velocities
#ax1.tricontourf(x, y, z, levels=40,linewidths=0.5, colors='k')
#cntr2 = ax1.tricontourf(x, y, z, levels=40, cmap="RdBu_r")
cm = plt.cm.get_cmap('RdYlBu')
#fig1.colorbar(cntr2, ax=ax1)

sc = ax2.scatter(conditional_array[~indices,1], conditional_array[~indices,2], marker='s',  s=4,c=velocities, cmap = cm)
#ax1.plot(x, max(y))
plt.colorbar(sc)
#ax1.scatter(conditional_array[indices,1], conditional_array[indices,2], c='red')
full_matrix = np.array([x,y,list_multipliers])
ax2.set_ylabel('Sweep offset [$\deg$]')
ax2.set_xlabel('Tail opening [$\deg$]')
ax2.grid(False)
plt.show()
fig2 = plt.figure()
ax2 = fig2.gca()
ax2.set_aspect('equal')
ax2.set_xlabel('$Re$', fontsize=20)  
ax2.set_ylabel('$Im$', fontsize=20)  
ax2.xaxis.set_tick_params(labelsize=10)
ax2.set_xlim(-1.1, 2.2)
ax2.set_ylim(-1.2, 1.2)
ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1.))
ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1.))
#ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.tick_params(labelsize=20)
circle = np.linspace(0,2*np.pi,101)
ax2.plot(np.cos(circle),np.sin(circle), linewidth=3.)
plt.grid(False)
ax2.scatter(full_matrix[2][0].real, full_matrix[2][0].imag , s=70,marker='.',color='red',
            facecolors='red', edgecolors='red', linewidth=1.8)
plt.gcf().subplots_adjust(left=0.16)