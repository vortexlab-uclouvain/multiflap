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

        damp = 1.0
        x0 = self.ms_obj.get_initial_guess()
        x = 1*(x0)
        ptlist = [x0]
        error = []
        dim = self.ms_obj.dim
        # Start the iterative scheme
        for i in range(self.max_iterations):

            start_timing = time.time()
            print("... Iteration i = " + str(i))

            MS, E, complete_solution = self.ms_obj.get_ms_scheme(x)

            epsilon = max(np.abs(E))
            error.append(epsilon)
            end_timing = time.time()
            print("... Iteration " + str(i) + " time: " + str(end_timing - start_timing))
            print("... Iteration number i = " + str(i) + " has error: " + str(norm(E, inf)))

    #        np.save(data_sim_history+"/Checkpoint_Iteration_"+str(i), x)
    #        np.save(data_sim_history+"/Error_History", error)


            if epsilon < self.tolerance:
                print("Iteration terminates at error : " + str(norm(E, inf)))
                print("Calculating Floquet Multipliers")

                ptlist.append(1*x)
                # jac_sg: jacobian semigroup (from one point to the next one)
                jac_sg = np.eye(dim)

                # Evaluation of Jacobian, using Semigroup (transition) property
                for i in range(0, self.ms_obj.point_number-1):
                    jac_sg = (MS[(i*dim):dim+(i*dim),
                                (i*dim):(i*dim)+dim]).dot(jac_sg)

                break

            else:

                # setting up LM equation 
                # [(MS.T)MS + lambda*diag((MS.T)MS)]*delta_x = (MS.T)E
                # put the Eq. in the form A*delta_x = B 

                A_0 = np.dot(MS.T,  MS)
                A = A_0 + damp* np.diag(np.diag(A_0))
                B = np.dot(MS.T, E)
                delta_x = np.linalg.solve(A, B)
                # Update of the solution

                for k in range(self.ms_obj.point_number):
                    x[k,:] = x[k,:] + delta_x[k*dim:(k+1)*dim]
                [E_new,_] = self.ms_obj.get_error_vector(x)

                if norm(E_new, inf) < norm(E, inf):
                    damp = damp / 10.0
                else :
                    damp = damp * 10.0
                ptlist.append(1*x)

        return x, ptlist, error, complete_solution, jac_sg
