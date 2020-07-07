import numpy as np
import bird_model as bm
from bird import Shoulder, Elbow, Wrist, Joint
from rk_integrator import rk2
from multiple_shooting import MultipleShooting
# construct the bird with BirdModel


bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                         axis_y=Joint(-np.deg2rad(19),np.pi/12,np.pi/2),
                         axis_z=Joint(0,np.deg2rad(42),np.pi))
bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                   axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                   axis_z=Joint(0.,0.,0.))

mybird = bm.BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist)

initial_time = 0.
final_time = 0.25
time_steps = 100
time_array = np.linspace(initial_time, final_time, time_steps)
x0 = [16, 0.4, 0.4, 0.4]
#lift = []
#for i in range (len(time_array)):
#    aeroforces = mybird.get_aeroforces(x0, time_array[i])
#    lift.append(aeroforces[1])
#
#Lift = np.asanyarray(lift)
#
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#
#fig = plt.figure()
#ax = fig.gca()
#plt.xlabel('1/T', fontsize=14)
#plt.ylabel('L(t)', fontsize=14)
#ax.plot(time_array, Lift, '-', label="Lift")
#plt.show()
#x0 = [16, 0, 0, 0]
#state_zero = mybird.dynamics(x0, 0)
point_distribution = MultipleShooting(x0, mybird).get_jacobian_analytical(x0, 0, 0.25)
