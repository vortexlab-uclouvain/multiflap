import numpy as np
import multiple_shooting as ms
import bird_model as bm
from bird import Shoulder, Elbow, Wrist, Joint
import time
from numpy.linalg import norm
from numpy import inf
import lma_solver as lma
bird_shoulder = Shoulder(axis_x=Joint(0.2,0.014,-np.pi/2),
                         axis_y=Joint(-np.deg2rad(19),np.deg2rad(20),np.pi/2),
                         axis_z=Joint(0,np.deg2rad(42),np.pi))
bird_elbow = Elbow(axis_x=Joint(0.,np.pi/6,-np.pi/2),
                   axis_y=Joint(np.pi/6,np.pi/6,-np.pi/2))

bird_wrist = Wrist(axis_y=Joint(-np.pi/6,np.pi/6,np.pi/2),
                   axis_z=Joint(0.,0.,0.))

mybird = bm.BirdModel(shoulder=bird_shoulder, elbow=bird_elbow, wrist=bird_wrist)
print('Bird obj created')
x0 = [16., 0.5, 0.1, 0.01]
ms_obj = ms.MultipleShooting(x0, mybird)
print('ms obj created')
mysol = lma.lma_solver(ms_obj = ms_obj).lma('/Users/gducci/UCL/PROJECT/Simulations/class_test')
