import numpy as np
import os
class SaveData:

    def __init__(self, simulation_name, folder_results = '../results'):
        self.simulation_name = simulation_name
        self.folder_results = folder_results
        self.folder_simulation =folder_results+'/'+simulation_name
    def save_data(self, filename, data_array):

        if filename == 'multipliers_LC':
                np.save(self.folder_simulation+'/'+'floquet_multipliers', data_array)
        elif filename == 'solution_LC':
                np.save(self.folder_simulation+'/'+'solution_LC', data_array)
        elif filename == 'time_array_LC':
                np.save(self.folder_simulation+'/'+'time_array', data_array)
        elif filename == 'fixed_points':
                np.save(self.folder_simulation+'/'+'fixed_points', data_array)
        else:
                np.save(self.folder_simulation+'/'+filename, data_array, allow_pickle=True)
        return print(filename+" file saved")

    def make_folder(self):
        path = self.folder_results+'/'+self.simulation_name
        try:
                os.mkdir(path)
                print("Folder correctly created")
        except:
                print("Error in trying to create data folder")

        return

    def remove_contents(self):
        path = self.folder_results+'/'+self.simulation_name
        try:
                os.remove(path+'/*')
        except:
                print("error removing file")

        return print("simulation folder created")


    def remove_directory(self):
        path = self.folder_results+'/'+self.simulation_name
        try:
                os.rmdir(path)
        except:
                print("error removing file")

        return ("Empty folder correctly removed")

    def write_simulation_settings(self, bird_obj):
        simulation_info = open(self.folder_simulation+'/simulation_settings.txt', 'w+')
        simulation_info.write('FRAMES COORDINATES (with respect to the aerodynamic reference system)'+'\n')
        simulation_info.write('Wing frame position: '+str(bird_obj.wingframe_position)+'\n')
        simulation_info.write('Tail frame position: '+str(bird_obj.wingframe_position_tail)+'\n')
        simulation_info.write('Tail opening: '+str(round(np.rad2deg(bird_obj.tail_opening)))+'\n')
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

    def write_simulation_results(self, error, multipliers):
        simulation_info = open(self.folder_simulation+'/DataFile_Results.txt', 'w+')
        simulation_info.write('Number of Iterations: '+str(np.size(error))+'\n')
        simulation_info.write('Error: '+str((error[-1]))+'\n')
        simulation_info.write('Floquet Multipliers: '+str(multipliers)+'\n')
        simulation_info.write('|$Lambda_1$|: '+str(np.linalg.norm(multipliers[0]))+'\n')
        simulation_info.write('|$Lambda_2$|: '+str(np.linalg.norm(multipliers[1]))+'\n')
        simulation_info.write('|$Lambda_3$|: '+str(np.linalg.norm(multipliers[2]))+'\n')
        simulation_info.write('|$Lambda_4$|: '+str(np.linalg.norm(multipliers[3]))+'\n')
        simulation_info.close() 
        
        return
