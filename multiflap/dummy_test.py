from utilities.LimitCycleForces import RetrievingAero
from aero_package.bird_model import BirdModel
import numpy as np
from aero_package.bird import Shoulder, Elbow, Wrist, Joint
from aero_package.settings import SimulationsSettings

bird_settings = SimulationsSettings(0.)
# generate bird kinematics by calling "bird" module
bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                       axis_y=Joint(-np.deg2rad(19),np.deg2rad(20),np.pi/2),
                       axis_z=Joint(0,np.deg2rad(42),np.pi))
bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                 axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                 axis_z=Joint(0.,0.,0.))

mybird = BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist, settings=bird_settings)

time_array = np.linspace(0, 0.25, 25)
u = np.ones(len(time_array))
w = np.zeros(len(time_array))
q = np.zeros(len(time_array))
theta = np.zeros(len(time_array))

x = np.array([u, w, q, theta]).T

retrieving = RetrievingAero(mybird, time_array)
