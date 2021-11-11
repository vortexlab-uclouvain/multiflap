from utilities.LimitCycleForces import RetrievingAero
from aero_package.bird_model import BirdModel
import numpy as np
from aero_package.bird import Shoulder, Elbow, Wrist, Joint
from aero_package.settings import SimulationsSettings
import matplotlib.pyplot as plt
amp_z = np.deg2rad(35)
bird_settings = SimulationsSettings(40.)
folder_name = '/Users/gducci/UCL/PROJECT/code/results/rerun/tail_45.000_sweep23.500'
# generate bird kinematics by calling "bird" module
bird_shoulder = Shoulder(axis_x=Joint(0.15, 0.0 ,-np.pi/2),
                         axis_y=Joint(-np.deg2rad(19),np.deg2rad(20),np.pi/2),
                         axis_z=Joint(0,np.deg2rad(42),np.pi))
bird_elbow = Elbow(axis_x=Joint(0.05,np.pi/6,-np.pi/2),
                   axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                   axis_z=Joint(0.,0.,0.))

mybird = BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist, settings=bird_settings)

time_array = np.linspace(0, 0.25, 50)
u = np.ones(len(time_array))*15
w = np.zeros(len(time_array))
q = np.zeros(len(time_array))
theta = np.zeros(len(time_array))

x = np.array([u, w, q, theta]).T
solution_LC = np.load(folder_name+ '/solution_LC.npy')
time_array_LC = np.load(folder_name+ '/time_array.npy')

retrieving = RetrievingAero(mybird, time_array)
aero = retrieving.postrun_aeroforces(x)
plt.plot(time_array, aero[1])
plt.plot(time_array, aero[2])
plt.show()
