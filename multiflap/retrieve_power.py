import numpy as np
from aero_package.bird_model import BirdModel
from aero_package.settings import SimulationsSettings
from aero_package.bird import Shoulder, Elbow, Wrist, Joint
from utilities.LimitCycleForces import RetrievingAero
import glob
import re
import os

def list_simulations(pathname):
    """
    pathname = /mypath/element_name
    example: pathname = './tail_*' to list all the elements starting for "tail_"
    """
    simlist = glob.glob(pathname)
    return simlist

sim = list_simulations('../results/power_check/tail_*_sweep_*_SAz_*')

def get_power(sim):
    
    list_results = []
    list_interpolation = []
    # String manipulation loop
    for i in range(len(sim)):
        simname = sim[i]
        
        # check wheter the limit cycle was found (this is rough!)
        # should put also a filter on the vertical velocity
        if os.path.isfile(simname+ '/solution_LC.npy'):
            x = re.findall(r'-?\d+\.?\d*', simname)
            x = list(map(float, x))
            
            solution = np.load(simname+'/solution_LC.npy' )
            sol_time = np.load(simname+'/time_array.npy')
            off_1 = np.load(simname+'/offset_elbow.npy')

            # Input of the simulation
            tail_opening = x[0]
            sweep_offset = x[1]
            amplitude = x[-1]
            # offset elbow
            f=open(simname+'/simulation_settings.txt')
            lines=f.readlines()
            line = lines[28]
            #print(line)
            off = re.findall(r'-?\d+\.?\d*', line)
            off = list(map(float, off))
            off = np.asanyarray(off)
            off = np.deg2rad(float(off))

            bird_settings = SimulationsSettings(tail_opening)
            bird_shoulder = Shoulder(axis_x=Joint(0.2, 0.014, -np.pi/2),
                                     axis_y=Joint(-np.deg2rad(sweep_offset),np.deg2rad(14),np.pi/2),
                                     axis_z=Joint(0,np.deg2rad(amplitude),np.pi))
            
            bird_elbow = Elbow(axis_x=Joint(off,np.pi/6,-np.pi/2),
                               axis_y=Joint(np.deg2rad(10),np.deg2rad(10),-np.pi/2))
            
            bird_wrist = Wrist(axis_y=Joint(-np.deg2rad(30),np.deg2rad(30),np.pi/2),
                               axis_z=Joint(0.,0.,0.))
            
            
            mybird = BirdModel(shoulder=bird_shoulder,
                               elbow=bird_elbow,
                               wrist=bird_wrist,
                               settings=bird_settings)
            
            post_running = RetrievingAero(mybird, sol_time)
            [_,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            P_ind,_,_] =  post_running.postrun_aeroforces(solution)
            
            np.save(simname+'/power_time_2.npy', P_ind)

write_power_2 = get_power(sim)
