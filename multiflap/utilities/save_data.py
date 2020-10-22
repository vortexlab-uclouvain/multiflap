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
        except:
                print("error")
        
        return print("simulation folder created")