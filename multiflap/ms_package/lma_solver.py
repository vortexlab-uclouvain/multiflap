import numpy as np
import time
from numpy.linalg import norm
from numpy import inf

class Solver:


    def __init__(self, tolerance = 1e-5, max_iterations = 100, ms_obj = None):
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.ms_obj = ms_obj

    def lma(self, result_directory):

        # Step zero of the scheme, storage of guessed points given as an input
        # by the user
        damp = 1.0
        #x_new = np.copy(x)
    #    data_sim_history = (result_directory+"/Simulation_History")
    #    os.mkdir(data_sim_history)
        x0 = self.ms_obj.get_initial_guess()
        x = 1*(x0)
        ptlist = [x0]
        error = []
        dim = self.ms_obj.dim
        # Start the iterative scheme
        for i in range(self.max_iterations):

            start_timing = time.time()
            print("... Iteration i = " + str(i))

            MS, dF, complete_solution = self.ms_obj.get_ms_scheme(x)

            epsilon = max(np.abs(dF))
            error.append(epsilon)
            end_timing = time.time()
            print("... Iteration " + str(i) + " time: " + str(end_timing - start_timing))
            print("... Iteration number i = " + str(i) + " has error: " + str(norm(dF, inf)))

    #        np.save(data_sim_history+"/Checkpoint_Iteration_"+str(i), x)
    #        np.save(data_sim_history+"/Error_History", error)


            if epsilon < self.tolerance:
                print("Iteration terminates at error : " + str(norm(dF, inf)))
                print("Calculating Floquet Multipliers")

                ptlist.append(1*x)
                Jacobian_semigroup = np.eye(dim)
                # =============================================================================
                #    Evaluation of Jacobian, using Semigroup property
                # =============================================================================
                for i in range(0, self.ms_obj.point_number- 1):
                    Jacobian_semigroup = (MS[(i*dim):dim+(i*dim),
                                             (i*dim):(i*dim)+dim]).dot(Jacobian_semigroup)

                break

            else:
                JJ = np.dot(MS.T,  MS)
                JdF = np.dot(MS.T, dF)
                H = JJ + damp* np.diag(np.diag(JJ))
                delta_x = np.linalg.solve(H, JdF)
                # Update of the solution

                for k in range(self.ms_obj.point_number):
                    x[k,:] = x[k,:] + delta_x[k*dim:(k+1)*dim]
                [dF_new,_] = self.ms_obj.get_error_vector(x)

                if norm(dF_new, inf) < norm(dF, inf):
                    damp = damp / 10.0
                else :
                    damp = damp * 10.0
                ptlist.append(1*x)

        return x, ptlist, error, complete_solution, Jacobian_semigroup
