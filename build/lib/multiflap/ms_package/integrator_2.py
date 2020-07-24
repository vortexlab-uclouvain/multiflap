import numpy as np
import collections


class Integrator:

    def __init__(self, ode_system, x0, time_array):

        self.ode_system = ode_system
        self.x0 = x0
        self.time_array = time_array
    def rk4(self):
        """
        Runge-Kutta 4 Integrator.
        Inputs:
        VelocityFunction: Function name to integrate
                          this function must have two inputs namely state space
                          vector and time. For example: velocity(ssp, t)
        InitialCondition: Initial condition, 1xd NumPy array, where d is the
                          dimension of the state space
        TimeArray: 1 x Nt NumPy array which contains instances for the solution
                   to be returned.
        Outputs:
        solution: Nt x d NumPy array which contains numerical solution of the
                       ODE.
        """
        rk_sol = collections.namedtuple('rk_sol',['x', 't'])
        #Generate the solution array to fill in:
        solution = np.zeros((np.size(self.time_array, 0),
                                  np.size(self.x0, 0)))
        #Assign the initial condition to the first element:
        solution[0, :] = self.x0

        for i in range(0, np.size(self.time_array) - 1):
            #Read time element:
            deltat = self.time_array[i + 1] - self.time_array[i]
            #Runge Kutta k's:
            k1 = deltat * self.ode_system(solution[i], self.time_array[i])
            k2 = deltat * self.ode_system(solution[i]+k1/2.0, self.time_array[i]+deltat/2.0)
            k3 = deltat * self.ode_system(solution[i]+k2/2.0, self.time_array[i]+deltat/2.0)
            k4 = deltat * self.ode_system(solution[i]+k3, self.time_array[i]+deltat)
            #Next integration step:
            solution[i + 1] = solution[i] + ((k1 +2*k2 + 2*k3 + k4)/6.0)

        sol = rk_sol(solution, self.time_array)
        return sol

    def rk3(self):
        """
        Runge-Kutta 3 Integrator.
        Inputs:
        VelocityFunction: Function name to integrate
                          this function must have two inputs namely state space
                          vector and time. For example: velocity(ssp, t)
        InitialCondition: Initial condition, 1xd NumPy array, where d is the
                          dimension of the state space
        TimeArray: 1 x Nt NumPy array which contains instances for the solution
                   to be returned.
        Outputs:
        solution: Nt x d NumPy array which contains numerical solution of the
                       ODE.
        """
        rk_sol = collections.namedtuple('rk_sol',['x', 't'])
        #Generate the solution array to fill in:
        solution = np.zeros((np.size(self.time_array, 0),
                                  np.size(self.x0, 0)))
        #Assign the initial condition to the first element:
        solution[0, :] = self.x0

        for i in range(0, np.size(self.time_array) - 1):
            #Read time element:
            deltat = self.time_array[i + 1] - self.time_array[i]

            #Runge Kutta k's:
            k1 = deltat * self.ode_system(solution[i], self.time_array[i])
            k2 = deltat * self.ode_system(solution[i]+k1/2.0, self.time_array[i]+deltat/2.0)
            k3 = deltat * self.ode_system(solution[i] -k1 + 2*k2, self.time_array[i]+deltat)
            #Next integration step:
            solution[i + 1] = solution[i] + ((k1 +4*k2 + k3)/6.0)

        sol = rk_sol(solution, self.time_array)
        return sol

    def rk2(self):
        """
        Runge-Kutta 2 Integrator.
        Inputs:
        VelocityFunction: Function name to integrate
                          this function must have two inputs namely state space
                          vector and time. For example: velocity(ssp, t)
        InitialCondition: Initial condition, 1xd NumPy array, where d is the
                          dimension of the state space
        TimeArray: 1 x Nt NumPy array which contains instances for the solution
                   to be returned.
        Outputs:
        solution: Nt x d NumPy array which contains numerical solution of the
                       ODE.
        """

        rk_sol = collections.namedtuple('rk_sol',['x', 't'])
        #Generate the solution array to fill in:
        solution = np.zeros((np.size(self.time_array, 0),
                                  np.size(self.x0, 0)))
        #Assign the initial condition to the first element:
        solution[0, :] = self.x0

        for i in range(0, np.size(self.time_array)-1):
            #Read time element:
            deltat = self.time_array[i + 1] - self.time_array[i]
            #Runge Kutta k's:
            k1 = deltat * self.ode_system(solution[i], self.time_array[i])
            k2 = deltat * self.ode_system(solution[i]+k1*(2./3.), self.time_array[i]+deltat*(2./3.))
            #Next integration step:
            solution[i + 1] = solution[i] + (k1/4. + (3./4.)*k2)

        sol = rk_sol(solution, self.time_array)
        return sol
