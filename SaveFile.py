import numpy as np

def SaveData(results_directory,ptlist, complete_solution,
             states_stack,
             error,
             eigenValues_SG,
             eigenVectors_SG,
             eigenValues,
             eigenVectors,
             Fx, Fy, Fz, Moment_total, F_tail, Moment_wing, Moment_tail, Moment_drag, Moment_lift):
    
    np.save(results_directory+'/outfile_ptlist', ptlist)
    np.save(results_directory+'/complete_solution', complete_solution)
    np.save(results_directory+'/final_points', complete_solution)
    np.save(results_directory+'/initial_points', states_stack)
    np.save(results_directory+'/Error', error)
    np.save(results_directory+'/outfile_JacobianEigenvalues_SG', eigenValues_SG)
    np.save(results_directory+'/outfile_JacobianEigenvector_SG', eigenVectors_SG)
    np.save(results_directory+'/outfile_JacobianEigenvalues', eigenValues)
    np.save(results_directory+'/outfile_JacobianEigenvector', eigenVectors)
    np.save(results_directory+'/Lift_coupled_v2', Fy)
    np.save(results_directory+'/Drag_coupled_v2', Fz)
    np.save(results_directory+'/Force_tail', F_tail)
    np.save(results_directory+'/Moment_total', Moment_total)
    np.save(results_directory+'/Moment_wing', Moment_wing)
    np.save(results_directory+'/Moment_lift', Moment_lift)
    np.save(results_directory+'/Moment_drag', Moment_drag)
    np.save(results_directory+'/Moment_tail', Moment_tail)
    



    
