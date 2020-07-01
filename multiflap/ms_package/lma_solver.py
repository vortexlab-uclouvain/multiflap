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
        lam= 1.0
        #x_new = np.copy(x)
    #    data_sim_history = (result_directory+"/Simulation_History")
    #    os.mkdir(data_sim_history)
        x0 = self.ms_obj.get_initial_guess()
        x = 1*(x0)
        ptlist = [x0]
        error = []
        # Start the iterative scheme
        for i in range(self.max_iterations):

            # For the step "i", evaluation of DF matrix and dF vector, calling the
            # MultiShootingMatrix function
            start_timing = time.time()
            print("... Iteration i = " + str(i))

            DF, dF, complete_solution = self.ms_obj.get_ms_scheme(x)

            # Evaluation of the maximum error in the dF error vector
            epsilon = max(np.abs(dF))
            error.append(epsilon)
            end_timing = time.time()
            print("... Iteration " + str(i) + " time: " + str(end_timing - start_timing))
            print("... Iteration number i = " + str(i) + " has error: " + str(norm(dF, inf)))

            # -------------------------------------------------------------------------
            # Save Data every iteration
            # -------------------------------------------------------------------------

    #        np.save(data_sim_history+"/Checkpoint_Iteration_"+str(i), x)
    #        np.save(data_sim_history+"/Error_History", error)


            # Refinement of the solution, by calculatin the dx vector
            if epsilon < self.tolerance:
                print("Iteration terminates at error : " + str(norm(dF, inf)))
                print("Calculating Floquet Multipliers")
    #            Jacobian = func.JacobianNumerical(x, 0, (M - 1)*tau)
                ptlist.append(1*x)
                Jacobian_semigroup = np.eye(self.ms_obj.dim)
                # =============================================================================
                #    Evaluation of Jacobian, using Semigroup property
                # =============================================================================
                for i in range(0, self.ms_obj.point_number- 1):
                    Jacobian_semigroup = (DF[(i*self.ms_obj.dim):self.ms_obj.dim+(i*self.ms_obj.dim),
                                             (i*self.ms_obj.dim):(i*self.ms_obj.dim)+self.ms_obj.dim]).dot(Jacobian_semigroup)

                break

            else:
                JJ = np.dot(DF.T,  DF)
                JdF = np.dot(DF.T, dF)
                H = JJ + lam* np.diag(np.diag(JJ))
                delta_x = np.linalg.solve(H, JdF)
                #delta_x = np.linalg.solve((DF), dF) # Matrice dx che trovo invertendo il sistema

                # Update of the solution
                for k in range(self.ms_obj.point_number):
                    x[k,:] = x[k,:] + delta_x[k*self.ms_obj.dim:(k+1)*self.ms_obj.dim]
                [dF_new,_] = self.ms_obj.get_df_vector(x)
                if norm(dF_new, inf) < norm(dF, inf):
                    #x = x_new;
                    # tau = tau
                    lam= lam / 10.0

                else :
                    lam= lam * 10.0
                    if lam> 1e10:
                        print ("lam is too large")
                ptlist.append(1*x)

        return x, ptlist, error, complete_solution, Jacobian_semigroup
