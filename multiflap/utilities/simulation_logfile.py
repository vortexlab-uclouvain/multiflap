import settings as settings
import numpy as np

'''
file simulation_write_data.py
@author Gianmarco Ducci
@copyright Copyright © UCLouvain 2020

multiflap is a Python tool for finding periodic orbits and assess their stability via the Floquet multipliers.

Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>

List of the contributors to the development of multiflap, Description and complete License: see LICENSE and NOTICE files.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
def write_simulation_settings(simulation_directory, bird_obj):
    simulation_info = open(simulation_directory+'/simulation_settings.txt', 'w+')
    simulation_info.write('FRAMES COORDINATES (with respect to the aerodynamic reference system)'+'\n')
    simulation_info.write('Wing frame position: '+str(bird_obj.wingframe_position)+'\n')
    simulation_info.write('Tail frame position: '+str(bird_obj.wingframe_position_tail)+'\n')
    simulation_info.write('Tail opening: '+str(round(np.rad2deg(bird_obj.tail_opening)))+'\n')
    simulation_info.write('Tail chord: '+str(settings.tail_length)+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('KINEMATICS'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('- SHOULDER JOINT:'+'\n')
    simulation_info.write('Amplitude shoulder x: '+str(np.rad2deg(bird_obj.shoulder.axis_x.amplitude))+'\n')
    simulation_info.write('Offset shoulder x: '+str(np.rad2deg(bird_obj.shoulder.axis_x.offset))+'\n')
    simulation_info.write('Phase shoulder x: '+str(np.rad2deg(bird_obj.shoulder.axis_x.phase))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude shoulder y: '+str(np.rad2deg(bird_obj.shoulder.axis_y.amplitude))+'\n')
    simulation_info.write('Offset shoulder y: '+str(np.rad2deg(bird_obj.shoulder.axis_y.offset))+'\n')
    simulation_info.write('Phase shoulder y: '+str(np.rad2deg(bird_obj.shoulder.axis_y.phase))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude shoulder z: '+str(np.rad2deg(bird_obj.shoulder.axis_z.amplitude))+'\n')
    simulation_info.write('Offset shoulder z: '+str(np.rad2deg(bird_obj.shoulder.axis_z.offset))+'\n')
    simulation_info.write('Phase shoulder z: '+str(np.rad2deg(bird_obj.shoulder.axis_z.phase))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('- ELBOW JOINT:'+'\n')
    simulation_info.write('Amplitude elbow y: '+str(np.rad2deg(bird_obj.elbow.axis_y.amplitude))+'\n')
    simulation_info.write('Offset elbow y: '+str(np.rad2deg(bird_obj.elbow.axis_y.offset))+'\n')
    simulation_info.write('Phase elbow y: '+str(np.rad2deg(bird_obj.elbow.axis_y.phase))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude elbow x: '+str(np.rad2deg(bird_obj.elbow.axis_x.amplitude))+'\n')
    simulation_info.write('Offset elbow x: '+str(np.rad2deg(bird_obj.elbow.axis_x.offset))+'\n')
    simulation_info.write('Phase elbow x: '+str(np.rad2deg(bird_obj.elbow.axis_x.phase))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('- WRIST JOINT:'+'\n')
    simulation_info.write('Amplitude wrist y: '+str(np.rad2deg(bird_obj.wrist.axis_y.amplitude))+'\n')
    simulation_info.write('Offset wrist y: '+str(np.rad2deg(bird_obj.wrist.axis_y.amplitude))+'\n')
    simulation_info.write('Phase wrist y: '+str(np.rad2deg(bird_obj.wrist.axis_y.amplitude))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude wrist z: '+str(np.rad2deg(bird_obj.wrist.axis_z.amplitude))+'\n')
    simulation_info.write('Offset wrist z: '+str(np.rad2deg(bird_obj.wrist.axis_z.amplitude))+'\n')
    simulation_info.write('Phase wrist z: '+str(np.rad2deg(bird_obj.wrist.axis_z.amplitude))+'\n')
    simulation_info.close() 

    return

def write_simulation_results(simulation_directory, error, multipliers):
    simulation_info = open(simulation_directory+'/DataFile_Results.txt', 'w+')
    simulation_info.write('Number of Iterations: '+str(np.size(error))+'\n')
    simulation_info.write('Error: '+str((error[-1]))+'\n')
    simulation_info.write('Floquet Multipliers: '+str(multipliers)+'\n')
    simulation_info.write('|$Lambda_1$|: '+str(np.linalg.norm(multipliers[0]))+'\n')
    simulation_info.write('|$Lambda_2$|: '+str(np.linalg.norm(multipliers[1]))+'\n')
    simulation_info.write('|$Lambda_3$|: '+str(np.linalg.norm(multipliers[2]))+'\n')
    simulation_info.write('|$Lambda_4$|: '+str(np.linalg.norm(multipliers[3]))+'\n')
    simulation_info.close() 
    
    return
