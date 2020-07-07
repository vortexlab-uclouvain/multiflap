import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from mpl_toolkits.mplot3d import proj3d
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import interp1d
import numpy.polynomial.polynomial as poly
import math
import LimitCycleForces as forces
import csv
import pandas as pd
"""
Nomenclature: To distinguish the coupled case from the non-coupled one we adopt
the following.
VariableName_wt = uncoupled case, wind tunnel (wt) configuration
VariableName_c  = coupled case
"""
#case_name_2 = 'RefinementTestCase26'
case_name = 'TestCase17'
with open('/Users/gducci/UCL/Simulations/utilities/Dataset2.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))
def readcsv(filename):
    data = pd.read_csv(filename) #Please add four spaces here before this line
    return(np.array(data))
yourArray = readcsv('/Users/gducci/UCL/Simulations/utilities/VictorData.csv')
yourArray[-1,1] = yourArray[0,1]
yourArray[-1,0] = 1.
retrieve_forces = False
results_directory = '/Users/gducci/UCL/Simulations/'+ case_name+'/Results'
if retrieve_forces == True:
    [Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift, _] = forces.ForceRetrieving_CaseName(case_name)
    np.save(results_directory+'/Lift_coupled_v2', Fy)
    np.save(results_directory+'/Drag_coupled_v2', Fz)
    np.save(results_directory+'/Force_tail', F_tail)
    np.save(results_directory+'/Moment_total', Moment_total)
    np.save(results_directory+'/Moment_wing', Moment_wing)
    np.save(results_directory+'/Moment_lift', Moment_lift)
    np.save(results_directory+'/Moment_drag', Moment_drag)
    np.save(results_directory+'/Moment_tail', Moment_tail)

checkpoint_directory = '/Users/gducci/UCL/Simulations/'+ case_name+'/Simulation_History'
#results_directory_2 = '/Users/gducci/UCL/Simulations/'+ case_name_2+'/Results'
#checkpoint_directory_2 = '/Users/gducci/UCL/Simulations/'+ case_name_2+'/Simulation_History'
#Figure_path = '/Users/gducci/UCL/Papers/Paper1_JNonlinearScience/PaperLatex/figures/results'
Figure_path = '/Users/gducci/UCL/presentiations/RF_meeting_1812'
# -------------------------------------------------------------------------
#  List of the filenames 
# -------------------------------------------------------------------------

jacobian_eigenvalues_filename = results_directory+'/outfile_JacobianEigenvalues.npy'
#jacobian_eigenvalues_filename_2 = results_directory_2+'/outfile_JacobianEigenvalues.npy'
periodic_orbit_filename = results_directory+'/complete_solution.npy'
list_point_filename = results_directory+'/outfile_ptlist.npy' 
error_filename = checkpoint_directory+'/Error_History.npy'
lift_filename = results_directory+'/Lift_coupled_v2.npy'
force_tail_filename = results_directory+'/Force_tail.npy'
drag_filename = results_directory+'/Drag_coupled_v2.npy'
Moment_total_filename = results_directory+'/Moment_total.npy'
Moment_wing_filename = results_directory+'/Moment_wing.npy'
Moment_tail_filename = results_directory+'/Moment_tail.npy'
Moment_drag_filename = results_directory+'/Moment_drag.npy'
Moment_lift_filename = results_directory+'/Moment_lift.npy'

jacobian_eigenvalues_filename_2 = results_directory+'/outfile_JacobianEigenvalues_SG.npy'
#periodic_orbit_filename_2 = results_directory_2+'/complete_solution.npy'
#list_point_filename_2 = results_directory_2+'/outfile_ptlist.npy' 
#error_filename_2 = checkpoint_directory_2+'/Error_History.npy'

# -------------------------------------------------------------------------
#  Open and import the file as an array
# -------------------------------------------------------------------------
#periodic_orbit = np.load(periodic_orbit_filename)
jacobian_eigenvalues = np.load(jacobian_eigenvalues_filename)   # Load jacobian eigenvalues   
jacobian_eigenvalues_analytical = np.load(jacobian_eigenvalues_filename_2)   # Load jacobian eigenvalues SemiGroup property                
error = np.load(error_filename)

#jacobian_eigenvalues_2 = np.load(jacobian_eigenvalues_filename_2)   # Load jacobian eigenvalues          
#error_2 = np.load(error_filename_2)
poly_grade = 16
periodic_orbit = np.load(periodic_orbit_filename)               # Load periodic orbit as an array
periodic_orbit = periodic_orbit.reshape(-1, periodic_orbit.shape[2])
lift_coupled = np.load(lift_filename)
drag_coupled = np.load(drag_filename)
force_tail = np.load(force_tail_filename)
Moment_total = np.load(Moment_total_filename)
Moment_wing = np.load(Moment_wing_filename)
Moment_tail = np.load(Moment_tail_filename)
Moment_drag = np.load(Moment_drag_filename)
Moment_lift = np.load(Moment_lift_filename)

time_array_SV = np.linspace(0, 1, len(lift_coupled))
velocity_U = periodic_orbit[:,0]
velocity_W = periodic_orbit[:,1]
Q = periodic_orbit[:,2]
Theta = periodic_orbit[:,3]
coef_lift = poly.polyfit(time_array_SV, lift_coupled, poly_grade)
coef_drag = poly.polyfit(time_array_SV, drag_coupled, poly_grade)
f_lift = poly.Polynomial(coef_lift)
f_drag = poly.Polynomial(coef_drag)

coeff_moment = poly.polyfit(time_array_SV, Moment_total, poly_grade)
coeff_moment_wing = poly.polyfit(time_array_SV, Moment_wing, poly_grade)

f_moment = poly.Polynomial(coeff_moment)
f_moment_wing = poly.Polynomial(coeff_moment_wing)
moment_fitted = f_moment(time_array_SV)
moment_fitted_wing = f_moment_wing(time_array_SV)
# =============================================================================
#  Flag quantities to plot
# =============================================================================
plot_eigenvalues = True
plot_error = True
plot_trajectory = True
plot_limit_cycle = True
plot_lift_BF = False
plot_drag_BF = False
plot_Fx_FF = True
plot_Fz_FF = True

# =============================================================================
#  Flag quantities to save
# =============================================================================
save_eigenvalues = False
save_error = False
save_trajectory = False
save_limit_cycle = False
save_lift_BF = False
save_drag_BF = False
save_Fx_FF = False
save_Fz_FF = False

# =============================================================================
# Plot The Lift. Wind tunnel vs. Coupled cases
# =============================================================================

#fig = plt.figure()
#ax = fig.gca() 
#fig.suptitle('Lift comparison', fontsize=18)
#plt.xlabel('1/T', fontsize=14)
#plt.ylabel('L(t)', fontsize=14)
#ax.plot(time_array_wt, lift_wt, '-.', label = "Non-coupled", linewidth=2.0)
#ax.plot(time_array_c, lift_c, '-.', color = 'red', label = "Coupled", linewidth=2.0)
#ax.legend()
#ax.set_xlim(0, 1)
#ax.grid(True)

# =============================================================================
# Plot eigenvalues
# =============================================================================
chord = 0.25
moment_norm = chord*1.2*9.81
fig1 = plt.figure()
ax1 = fig1.gca() 
#fig1.suptitle('Moments', fontsize=18)
plt.xlabel('1/T', fontsize=14)
plt.ylabel('$M(t)/(m\cdot g\cdot c)$', fontsize=14)
ax1.plot(time_array_SV, moment_fitted_wing/moment_norm, '-', label = "Moment wing", linewidth = 1., color = 'red')
ax1.plot(time_array_SV, Moment_tail/moment_norm, '-', label = "Moment tail", linewidth = 1., color = 'green')
ax1.plot(time_array_SV, moment_fitted/moment_norm, '-', label = "Total moment", linewidth = 2., color = 'blue')
#ax1.plot(yourArray[:,0], yourArray[:,1], '.-', label = "Wing only + control", linewidth = 1., color = 'red')

ax1.legend()
ax1.set_xlim(0, 1)
ax1.grid(True)
plt.savefig(Figure_path+'/Moment_me.eps', format = 'eps')

moment_fit = f_drag(time_array_SV)
fig1 = plt.figure()
ax1 = fig1.gca() 
fig1.suptitle('Moments', fontsize=18)
plt.xlabel('1/T', fontsize=14)
plt.ylabel('$M_i(t)$', fontsize=14)
ax1.plot(time_array_SV, Moment_drag, '-', label = "Moment drag")
ax1.plot(time_array_SV, Moment_lift, '-', label = "Moment lift")
ax1.plot(time_array_SV, Moment_tail, '-', label = "Moment tail")

ax1.legend()
ax1.set_xlim(0, 1)
ax1.grid(True)

if plot_eigenvalues==True:
    fig5 = plt.figure()
    ax5 = fig5.gca()
    ax5.set_aspect('equal')
    ax5.set_xlabel('$Re$', fontsize=20)  
    ax5.set_ylabel('$Im$', fontsize=20)  
    ax5.xaxis.set_tick_params(labelsize=10)
    ax5.set_xlim(-1.1, 1.1)
    ax5.set_ylim(-1.1, 1.1)
    ax5.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1.))
    ax5.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1.))
    ax5.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tick_params(labelsize=20)
    circle = np.linspace(0,2*np.pi,101)
    ax5.plot(np.cos(circle),np.sin(circle), linewidth=1., color='red')
#    plt.grid( linestyle='-', linewidth=0.4)
    plt.grid(True)
    ax5.scatter(jacobian_eigenvalues.real, jacobian_eigenvalues.imag , s=50,marker='.',facecolors='none', edgecolors='blue', linewidth=2.)
#    ax5.scatter(jacobian_eigenvalues_analytical.real, jacobian_eigenvalues_analytical.imag , s=20,marker='x',facecolors='black', edgecolors='black', linewidth=2.,label='M = 2')
    plt.gcf().subplots_adjust(left=0.1)
    plt.gcf().subplots_adjust(bottom=0.16)

    if save_eigenvalues == True:
        plt.savefig(Figure_path+'/jacobian_eigenvalues.eps', format = 'eps')

# =============================================================================
# Plot error
# =============================================================================

if plot_error == True:
    fig6 = plt.figure()
    ax6 = fig6.gca()
#    ax6.set_aspect('equal')
    ax6.set_yscale('log')
    ax6.set_xlabel('Iteration Number', fontsize=18)  
    ax6.set_ylabel('Error', fontsize=18)  
    ax6.xaxis.set_tick_params(labelsize=10)
    ax6.plot(error, '-o')
    #ax6.plot(error_2, '-o', label='$M=8$')
    ax6.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    plt.grid()
    if save_error == True:
        plt.savefig(Figure_path+'/error.eps', format = 'eps')

# =============================================================================
# Plot error
# =============================================================================

#if plot_error==True:
#    fig7 = plt.figure()
#    ax7 = fig7.gca()
##    ax7.set_aspect('equal')
#    ax7.set_yscale('log')
#    ax7.set_xlabel('Iteration Number', fontsize=18)  
#    ax7.set_ylabel('Error', fontsize=18)  
#    ax7.xaxis.set_tick_params(labelsize=10)
#    #ax7.plot(error_2, '-o')
#    ax7.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#    plt.grid()
#    if save_error == True:
#        plt.savefig(Figure_path+'/error.eps', format = 'eps')


# =============================================================================
# Plot error
# =============================================================================
        
u_hor = velocity_U*np.cos(Theta) + velocity_W*np.sin(Theta)
v_vert = velocity_U*np.sin(Theta) - velocity_W*np.cos(Theta)
z_initial = 10
x = ((np.cumsum(u_hor))/len(u_hor))/4
z = z_initial + (np.cumsum(v_vert))/len(u_hor)

if plot_trajectory==True:
    fig8 = plt.figure()  
    ax8 = fig8.gca() 
    plt.xlabel('X [m]', fontsize=14)
    plt.ylabel('Z [m]', fontsize=14)
#    ax8.set_ylim(9.8, 10)
    #ax9.set_xlim(0, max(x))
    ax8.grid(True)
    ax8.plot(x, z, color = 'green',linewidth=3)
    
    plt.gcf().subplots_adjust(left=0.15)
    plt.show()
    if save_trajectory == True:
        plt.savefig(Figure_path+'/X-Z_trajectory.eps', format = 'eps')

# =============================================================================
# Plot limit cycle
# =============================================================================
    
if plot_limit_cycle==True:
    def contour(k_min, k_max, Theta, ThetaMin, ThetaMax):
        color = ((k_min - k_max)/(ThetaMin - ThetaMax))*Theta + k_min - ((k_min - k_max)/(ThetaMin - ThetaMax))*ThetaMin
        return color
    
    N = len(Theta)
    
    spacing_U = np.abs((np.abs(max(velocity_U)) - np.abs(min(velocity_U))))/3
    spacing_W = np.abs((np.abs(max(velocity_W)) - np.abs(min(velocity_W))))/3
    spacing_Q = (np.abs(max(Q)) - np.abs(min(Q)))*5
    
    ThetaMax = (max(Theta))
    ThetaMin = (min(Theta))
    fig9 = plt.figure()
    ax9 = fig9.gca(projection='3d')
    ax9.view_init( 25, -135)
    ax9.w_xaxis.set_pane_color((1., 0.95, 0.95, 0.9))
    ax9.w_yaxis.set_pane_color((1., 0.95, 0.95, 0.9))
    ax9.w_zaxis.set_pane_color((1., 0.95, 0.95, 0.9))
    ax9.set_xlabel('$u$', fontsize=20)
    ax9.set_ylabel('$w$', fontsize=20)
    ax9.set_zlabel('$q$', fontsize=20)
    ax9.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax9.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax9.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.tick_params(labelsize=12)
#    ax9.xaxis.set_major_locator(mpl.ticker.MultipleLocator(spacing_U))
#    ax9.yaxis.set_major_locator(mpl.ticker.MultipleLocator(spacing_W))
#    ax9.zaxis.set_major_locator(mpl.ticker.MultipleLocator(spacing_Q))
    
    #sc = ax10.scatter(velocity_U,velocity_W, Q,
    #             c=mpl.cm.jet(contour(0.2, 0.9, Theta, ThetaMin, ThetaMax)), edgecolor='none')
    
    sc = ax9.scatter(velocity_U,velocity_W, Q,
                 c=Theta, cmap = 'jet', edgecolor='none')
    
    cb = fig9.colorbar(sc)
    cb.ax.tick_params(labelsize=15)
    cb.set_label('Theta Angle, $\Theta$', labelpad=-60, y=1.1, rotation=0, fontsize=18)
    if save_limit_cycle == True:
        plt.savefig(Figure_path+'/LimitCycle.eps', format = 'eps')



# =============================================================================
# Plot Forces in the body frame
# =============================================================================

if plot_lift_BF == True:
    lift_BF = f_lift(time_array_SV)
    time_array_wt = np.linspace(0, 1, len(lift_coupled))
    fig10 = plt.figure()
    ax10 = fig10.gca() 
    y_lim_inf = min(lift_BF)
    y_lim_sup = max(lift_BF)
    
    if y_lim_inf > 0:
        y_lim_inf = math.ceil((y_lim_inf))
    else:
        y_lim_inf = math.floor((y_lim_inf))
    
    if y_lim_sup > 0:
        y_lim_sup = math.ceil((y_lim_sup))
    else:
        y_lim_sup = math.floor((y_lim_sup))

    ax10.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.25))
    ax10.yaxis.set_major_locator(mpl.ticker.MultipleLocator(3.))
    ax10.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax10.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tick_params(labelsize=20)

    plt.xlabel('t/T', fontsize=20)
    plt.ylabel('$L(t)$', fontsize=20)
    ax10.plot(time_array_SV, lift_BF, linewidth=3.0)
#    ax11.plot(time_array_SV, lift_coupled, label = "Lift", linewidth=0.5, color='red')
    ax10.set_xlim(0, 1)
    ax10.grid(True)
    plt.gcf().subplots_adjust(left=0.165)
    plt.gcf().subplots_adjust(bottom=0.15)

    if save_lift_BF == True:
        plt.savefig(Figure_path+'/lift_BF.eps', format = 'eps')
    plt.show()

if plot_drag_BF == True:
    drag_BF = f_drag(time_array_SV)
    time_array_wt = np.linspace(0, 1, len(lift_coupled))
    fig11 = plt.figure()
    ax11 = fig11.gca() 
    y_lim_inf = min(drag_BF)
    y_lim_sup = max(drag_BF)
    
    if y_lim_inf > 0:
        y_lim_inf = math.ceil((y_lim_inf))
    else:
        y_lim_inf = math.floor((y_lim_inf))
    
    if y_lim_sup > 0:
        y_lim_sup = math.ceil((y_lim_sup))
    else:
        y_lim_sup = math.floor((y_lim_sup))

    ax11.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.25))
    ax11.yaxis.set_major_locator(mpl.ticker.MultipleLocator(2.))
    ax11.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax11.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tick_params(labelsize=20)

    plt.xlabel('t/T', fontsize=20)
    plt.ylabel('$D(t)$', fontsize=20)
    ax11.plot(time_array_SV, drag_BF, linewidth=3.0)
#    ax11.plot(time_array_SV, drag_coupled, label = "Drag", linewidth=0.5, color='red')
    ax11.set_xlim(0, 1)
    ax11.grid(True)
    plt.gcf().subplots_adjust(left=0.155)
    plt.gcf().subplots_adjust(bottom=0.15)

    if save_drag_BF == True:
        plt.savefig(Figure_path+'/drag_BF.eps', format = 'eps')
    plt.show()


# =============================================================================
# Plot Forces in the inertial frame
# =============================================================================
u_dot = -Q*velocity_W - 9.81*np.sin(Theta) - drag_coupled/1.2
w_dot = Q*velocity_U + 9.81*np.cos(Theta) - (lift_coupled + force_tail)/1.2
Fx = (u_dot*np.cos(Theta) + w_dot*np.sin(Theta))*1.2
Fz = (-(w_dot*np.cos(Theta) - u_dot*np.sin(Theta))*1.2)
print(np.mean(Fz))
coef_Fx_FF = poly.polyfit(time_array_SV, Fx, poly_grade)
coef_Fz_FF = poly.polyfit(time_array_SV, Fz, poly_grade)
f_Fx = poly.Polynomial(coef_Fx_FF)
f_Fz = poly.Polynomial(coef_Fz_FF)

Fx_FF = f_Fx(time_array_SV)
Fz_FF = f_Fz(time_array_SV)
force_mean_Fx = [np.mean(Fx_FF/(1.2*9.81))]*len(time_array_SV)
force_mean_Fz = [np.mean(Fz_FF)]*len(time_array_SV)

if plot_Fx_FF == True:
    
    
    fig12 = plt.figure()
    ax12 = fig12.gca()
    plt.xlabel('t/T', fontsize=20)
    plt.ylabel('$F_{hor}/w$', fontsize=20)
    y_lim_inf = min(Fx_FF)
    y_lim_sup = max(Fx_FF)
    
    if y_lim_inf > 0:
        y_lim_inf = math.ceil((y_lim_inf))
    else:
        y_lim_inf = math.floor((y_lim_inf))
    
    if y_lim_sup > 0:
        y_lim_sup = math.ceil((y_lim_sup))
    else:
        y_lim_sup = math.floor((y_lim_sup))

#    ax12.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.25))
#    ax12.yaxis.set_major_locator(mpl.ticker.MultipleLocator(3.))
    ax12.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax12.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tick_params(labelsize=20)
    ax12.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.25))
    ax12.plot(time_array_SV, force_mean_Fx ,linewidth=2., color='red', label = 'Mean value')
    ax12.plot(time_array_SV, Fx_FF/(1.2*9.81) ,linewidth=3., label = 'Mean value')
    ax12.set_xlim(0, 1)
    ax12.grid(True)
    plt.gcf().subplots_adjust(left=0.155)
    plt.gcf().subplots_adjust(bottom=0.15)

    if save_Fx_FF == True:
        plt.savefig(Figure_path+'/Fx_FF.eps', format = 'eps')
    plt.show()

if plot_Fz_FF == True:
    
    
    fig13 = plt.figure()
    ax13 = fig13.gca()
    plt.xlabel('t/T', fontsize=20)
    plt.ylabel('$F_{ver}/w$', fontsize=20)
#    y_lim_inf = min(Fz_FF)
#    y_lim_sup = max(Fz_FF)
#    
#    if y_lim_inf > 0:
#        y_lim_inf = math.ceil((y_lim_inf))
#    else:
#        y_lim_inf = math.floor((y_lim_inf))
#    
#    if y_lim_sup > 0:
#        y_lim_sup = math.ceil((y_lim_sup))
#    else:
#        y_lim_sup = math.floor((y_lim_sup))
    
#    ax13.set_ylim(y_lim_inf, y_lim_sup)
#    ax13.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.25))
#    ax13.yaxis.set_major_locator(mpl.ticker.MultipleLocator(3.))
    plt.tick_params(labelsize=16)
    ax13.xaxis.set_major_locator(mpl.ticker.MultipleLocator(.25))
    ax13.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax13.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

#    ax13.plot(time_array_SV, force_mean_Fx ,linewidth=2., color='red', label = 'Mean value')
    ax13.plot(time_array_SV, force_mean_Fz ,linewidth=2., color='red', label = 'Mean value')
    ax13.plot(time_array_SV, Fz_FF/(1.2*9.81) ,linewidth=3., label = 'Mean value')
    ax13.plot(time_array_SV, Fz/(1.2*9.81) ,linewidth=3., label = 'Mean value')
    ax13.set_xlim(0, 1)
    ax13.grid(True)
    plt.gcf().subplots_adjust(left=0.14)
    plt.gcf().subplots_adjust(bottom=0.15)

    if save_Fz_FF == True:
        plt.savefig(Figure_path+'/Fz_FF.eps', format = 'eps')
    plt.show()
