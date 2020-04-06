# -*- coding: utf-8 -*-
import settings as settings
import numpy as np

def DataFileParameters(simulation_directory):
    simulation_info = open(simulation_directory+'/DataFile_Parameters.txt', 'w+')
    simulation_info.write('Period: '+str(1/(settings.frequency))+'\n')
    simulation_info.write('FRAMES COORDINATES (with respect to the aerodynamic reference system)'+'\n')
    simulation_info.write('Wing frame position: '+str(settings.wingframe_position)+'\n')
    simulation_info.write('Tail frame position: '+str(settings.wingframe_position_tail)+'\n')
    simulation_info.write('Tail opening: '+str(settings.tail_opening)+'\n')
    simulation_info.write('Tail chord: '+str(settings.tail_length)+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('KINEMATICS'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('- SHOULDER JOINT:'+'\n')
    simulation_info.write('Amplitude shoulder x: '+str(np.rad2deg(settings.amplitude_shoulder_x))+'\n')
    simulation_info.write('Offset shoulder x: '+str(np.rad2deg(settings.offset_shoulder_x))+'\n')
    simulation_info.write('Phase shoulder x: '+str(np.rad2deg(settings.phaseangle_shoulder_x))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude shoulder y: '+str(np.rad2deg(settings.amplitude_shoulder_y))+'\n')
    simulation_info.write('Offset shoulder y: '+str(np.rad2deg(settings.offset_shoulder_y))+'\n')
    simulation_info.write('Phase shoulder y: '+str(np.rad2deg(settings.phaseangle_shoulder_y))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude shoulder z: '+str(np.rad2deg(settings.amplitude_shoulder_z))+'\n')
    simulation_info.write('Offset shoulder z: '+str(np.rad2deg(settings.offset_shoulder_z))+'\n')
    simulation_info.write('Phase shoulder z: '+str(np.rad2deg(settings.phaseangle_shoulder_z))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('- ELBOW JOINT:'+'\n')
    simulation_info.write('Amplitude elbow y: '+str(np.rad2deg(settings.amplitude_elbow_y))+'\n')
    simulation_info.write('Offset elbow y: '+str(np.rad2deg(settings.offset_elbow_y))+'\n')
    simulation_info.write('Phase elbow y: '+str(np.rad2deg(settings.phaseangle_elbow_y))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude elbow x: '+str(np.rad2deg(settings.amplitude_elbow_x))+'\n')
    simulation_info.write('Offset elbow x: '+str(np.rad2deg(settings.offset_elbow_x))+'\n')
    simulation_info.write('Phase elbow x: '+str(np.rad2deg(settings.phaseangle_elbow_x))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('- WRIST JOINT:'+'\n')
    simulation_info.write('Amplitude wrist y: '+str(np.rad2deg(settings.amplitude_wrist_y))+'\n')
    simulation_info.write('Offset wrist y: '+str(np.rad2deg(settings.offset_wrist_y))+'\n')
    simulation_info.write('Phase wrist y: '+str(np.rad2deg(settings.phaseangle_wrist_y))+'\n')
    simulation_info.write('----------------------------------'+'\n')
    simulation_info.write('Amplitude wrist z: '+str(np.rad2deg(settings.amplitude_wrist_z))+'\n')
    simulation_info.write('Offset wrist z: '+str(np.rad2deg(settings.offset_wrist_z))+'\n')
    simulation_info.write('Phase wrist z: '+str(np.rad2deg(settings.phaseangle_wrist_z))+'\n')
    simulation_info.close() 

    return

def DataFileResults(simulation_directory, error, multipliers):
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


    
    