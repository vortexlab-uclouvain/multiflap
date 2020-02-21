import numpy as np
import matplotlib.pyplot as plt

"""
Nomenclature: To distinguish the coupled case from the non-coupled one we adopt
the following.
VariableName_wt = uncoupled case, wind tunnel (wt) configuration
VariableName_c  = coupled case
"""
latex_image_folder = '/Users/gducci/UCL/Report/Meeting180619/figures/liftingline/'

# -------------------------------------------------------------------------
# List of the postprocessing filenames 
# -------------------------------------------------------------------------

lift_wt_filename = 'outfile_lift_u15.npy'       # Lift wind tunnel condition
lift_c_filename = 'outfile_lift_coupled_2.npy'  # Lift coupled with dynamics
velocity_U_filename = 'outfile_Uvelocity.npy'   # U-velocity state variable
velocity_W_filename = 'outfile_Wvelocity.npy'   # W-velocity state variable
Q_filenane = 'outfile_Q.npy'                    # Q-ang. velocity state var.
Theta_filename = 'outfile_Theta.npy'            # Theta angle state variable

# -------------------------------------------------------------------------
# Open and import the file as an array
# -------------------------------------------------------------------------

lift_wt = np.load(lift_wt_filename)             
lift_c = np.load(lift_c_filename)
velocity_U = np.load(velocity_U_filename)
velocity_W = np.load(velocity_W_filename)
Q = np.load(Q_filenane)
Theta = np.load(Theta_filename)

time_array_wt = np.linspace(0, 1, len(lift_wt))
time_array_c = np.linspace(0, 1, len(lift_c))
time_array_vel = np.linspace(0, 1, len(velocity_U))

# -------------------------------------------------------------------------
# ● Plot The Lift. Wind tunnel vs. Coupled cases
# -------------------------------------------------------------------------

fig = plt.figure()
ax = fig.gca() 
fig.suptitle(None, fontsize=18)
plt.xlabel('1/T', fontsize=14)
plt.ylabel('L(t)', fontsize=14)
ax.plot(time_array_wt, lift_wt,  label = "Non-coupled", linewidth=2.0)
ax.plot(time_array_c, lift_c, color = 'red', label = "Coupled", linewidth=2.0)
ax.legend()
ax.set_xlim(0, 1)
ax.grid(True)
plt.savefig(latex_image_folder + 'coupling_comparison.png', dpi = 600)

# -------------------------------------------------------------------------
# ● Plot U-velocity
# -------------------------------------------------------------------------

fig1 = plt.figure()
ax1 = fig1.gca() 
fig1.suptitle('U velocity', fontsize=18)
plt.xlabel('Time 1/T', fontsize=14)
plt.ylabel('U(t)', fontsize=14)
ax1.set_xlim(0, 1)
ax1.grid(True)
ax1.plot(time_array_vel,  velocity_U, '.', color = 'green')

# -------------------------------------------------------------------------
# ● Plot W-velocity
# -------------------------------------------------------------------------

fig2 = plt.figure()  
ax2 = fig2.gca() 
fig2.suptitle('W velocity', fontsize=18)
plt.xlabel('Time 1/T', fontsize=14)
plt.ylabel('W(t)', fontsize=14)
ax2.set_xlim(0, 1)
ax2.grid(True)
ax2.plot(time_array_vel,  velocity_W, '.', color = 'red')

# -------------------------------------------------------------------------
# ● Plot Q
# -------------------------------------------------------------------------

fig3 = plt.figure()  
ax3 = fig3.gca() 
fig3.suptitle('Q angular velocity', fontsize=18)
plt.xlabel('Time 1/T', fontsize=14)
plt.ylabel('Q(t)', fontsize=14)
ax3.set_xlim(0, 1)
ax3.grid(True)
ax3.plot(time_array_vel,  Q, '.', color = 'blue')

# -------------------------------------------------------------------------
# ● Plot Theta
# -------------------------------------------------------------------------

fig4 = plt.figure()  
ax4 = fig4.gca() 
fig4.suptitle('Theta angle', fontsize=18)
plt.xlabel('Time 1/T', fontsize=14)
plt.ylabel('$\Theta$(t)', fontsize=14)
ax4.set_xlim(0, 1)
ax4.grid(True)
ax4.plot(time_array_vel,  Theta, '.', color = 'red')

plt.show()