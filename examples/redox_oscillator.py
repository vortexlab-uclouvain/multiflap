import numpy as np
import multiflap as mf
import matplotlib.pyplot as plt

x = [0.5, 0.5, 0.6, 0.2]

time_array = np.linspace(0, 180, 90000)
mymodel = mf.RedoxModel()

ms_obj =  mf.MultipleShootingPeriod(x, M=2, period_guess= 23.,
                                    t_steps=90000, model=mymodel, integrator='odeint')

mysol = mf.SolverPeriod(ms_obj = ms_obj).lma()

jac = mysol[4]

eigenvalues, eigenvectors = np.linalg.eig(jac)

sol_array = mysol[3].space
sol_time = mysol[3].time
period = sol_time[-1]

# Save simulation data
sim_name = 'study_case'
mf.SaveData(sim_name, folder_results = './results').make_folder()
mf.SaveData(sim_name).save_data('state_space_solution', sol_array)

# Fast plot simulation results
plot = mf.Plot()
plot.limit_cycle_2D(sol_array[:,0], sol_array[:,1])
plot.limit_cycle_3D(sol_array[:,0], sol_array[:,1], sol_array[:,3])
plot.plot_multipliers(eigenvalues)
label = ['D1', 'D2', 'R', 'A']
plot.plot_time_series(sol_array, sol_time, label=label)
