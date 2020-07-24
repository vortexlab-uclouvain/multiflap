'''
file rk_integrator.py
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
import numpy as np
import collections

def rk4(ode_system, x0, time_array):
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
    solution = np.zeros((np.size(time_array, 0),
                              np.size(x0, 0)))
    #Assign the initial condition to the first element:
    solution[0, :] = x0

    for i in range(0, np.size(time_array) - 1):
        #Read time element:
        deltat = time_array[i + 1] - time_array[i]
        #Runge Kutta k's:
        k1 = deltat * ode_system(solution[i], time_array[i])
        k2 = deltat * ode_system(solution[i]+k1/2.0, time_array[i]+deltat/2.0)
        k3 = deltat * ode_system(solution[i]+k2/2.0, time_array[i]+deltat/2.0)
        k4 = deltat * ode_system(solution[i]+k3, time_array[i]+deltat)
        #Next integration step:
        solution[i + 1] = solution[i] + ((k1 +2*k2 + 2*k3 + k4)/6.0)

    sol = rk_sol(solution, time_array)
    return sol

def rk3(ode_system, x0, time_array):
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
    solution = np.zeros((np.size(time_array, 0),
                              np.size(x0, 0)))
    #Assign the initial condition to the first element:
    solution[0, :] = x0

    for i in range(0, np.size(time_array) - 1):
        #Read time element:
        deltat = time_array[i + 1] - time_array[i]

        #Runge Kutta k's:
        k1 = deltat * ode_system(solution[i], time_array[i])
        k2 = deltat * ode_system(solution[i]+k1/2.0, time_array[i]+deltat/2.0)
        k3 = deltat * ode_system(solution[i] -k1 + 2*k2, time_array[i]+deltat)
        #Next integration step:
        solution[i + 1] = solution[i] + ((k1 +4*k2 + k3)/6.0)

    sol = rk_sol(solution, time_array)
    return sol

def rk2( ode_system, x0, time_array):
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
    solution = np.zeros((np.size(time_array, 0),
                              np.size(x0, 0)))
    #Assign the initial condition to the first element:
    solution[0, :] = x0

    for i in range(0, np.size(time_array)-1):
        #Read time element:
        deltat = time_array[i + 1] - time_array[i]
        #Runge Kutta k's:
        k1 = deltat * ode_system(solution[i], time_array[i])
        k2 = deltat * ode_system(solution[i]+k1*(2./3.), time_array[i]+deltat*(2./3.))
        #Next integration step:
        solution[i + 1] = solution[i] + (k1/4. + (3./4.)*k2)

    sol = rk_sol(solution, time_array)
    return sol
