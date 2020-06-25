import numpy as np  # Import numpy
#import flapping_forces as FF
import scipy.integrate as ode
import bird_model as bm
def RK4(velocityFunction, initialCondition, timeArray, **optional):
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
    SolutionArray: Nt x d NumPy array which contains numerical solution of the
                   ODE.
    """
    #Generate the solution array to fill in:
    SolutionArray = np.zeros((np.size(timeArray, 0),
                              np.size(initialCondition, 0)))
    #Assign the initial condition to the first element:
    SolutionArray[0, :] = initialCondition

    for i in range(0, np.size(timeArray) - 1):
        #Read time element:
        deltat = timeArray[i + 1] - timeArray[i]
        #Runge Kutta k's:
        k1 = deltat * velocityFunction(SolutionArray[i], timeArray[i], **optional)
        k2 = deltat * velocityFunction(SolutionArray[i]+k1/2.0, timeArray[i]+deltat/2.0, **optional)
        k3 = deltat * velocityFunction(SolutionArray[i]+k2/2.0, timeArray[i]+deltat/2.0, **optional)
        k4 = deltat * velocityFunction(SolutionArray[i]+k3, timeArray[i]+deltat, **optional)
        #Next integration step:
        SolutionArray[i + 1] = SolutionArray[i] + ((k1 +2*k2 + 2*k3 + k4)/6.0)
    return SolutionArray

def RK3(velocityFunction, initialCondition, timeArray, **optional):
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
    SolutionArray: Nt x d NumPy array which contains numerical solution of the
                   ODE.
    """
    #Generate the solution array to fill in:
    SolutionArray = np.zeros((np.size(timeArray, 0),
                              np.size(initialCondition, 0)))
    #Assign the initial condition to the first element:
    SolutionArray[0, :] = initialCondition

    for i in range(0, np.size(timeArray) - 1):
        #Read time element:
        deltat = timeArray[i + 1] - timeArray[i]
        
        #Runge Kutta k's:
        k1 = deltat * velocityFunction(SolutionArray[i], timeArray[i], **optional)
        k2 = deltat * velocityFunction(SolutionArray[i]+k1/2.0, timeArray[i]+deltat/2.0, **optional)
        k3 = deltat * velocityFunction(SolutionArray[i] -k1 + 2*k2, timeArray[i]+deltat, **optional)
        #Next integration step:
        SolutionArray[i + 1] = SolutionArray[i] + ((k1 +4*k2 + k3)/6.0)
        
    return SolutionArray

def RK2( velocityFunction, initialCondition, timeArray, **optional):
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
    SolutionArray: Nt x d NumPy array which contains numerical solution of the
                   ODE.
    """
    #Generate the solution array to fill in:
    SolutionArray = np.zeros((np.size(timeArray, 0),
                              np.size(initialCondition, 0)))
    #Assign the initial condition to the first element:
    SolutionArray[0, :] = initialCondition

    for i in range(0, np.size(timeArray)-1):
        #Read time element:
        deltat = timeArray[i + 1] - timeArray[i]
        #Runge Kutta k's:
        k1 = deltat * velocityFunction(SolutionArray[i], timeArray[i], **optional)
        k2 = deltat * velocityFunction(SolutionArray[i]+k1*(2/3), timeArray[i]+deltat*(2./3.), **optional)
        #Next integration step:
        SolutionArray[i + 1] = SolutionArray[i] + (k1/4. + (3./4.)*k2)
        
    return SolutionArray

def RK2_Gust(velocityFunction, initialCondition, timeArray, **optional):
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
    SolutionArray: Nt x d NumPy array which contains numerical solution of the
                   ODE.
    """
    #Generate the solution array to fill in:
    SolutionArray = np.zeros((np.size(timeArray, 0),
                              np.size(initialCondition, 0)))
    #Assign the initial condition to the first element:
    SolutionArray[0, :] = initialCondition

    for i in range(0, np.size(timeArray)-1):
        #Read time element:
        deltat = timeArray[i + 1] - timeArray[i]
        #Runge Kutta k's:
        k1 = deltat * velocityFunction(SolutionArray[i], timeArray[i], **optional)[0]
        k2 = deltat * velocityFunction(SolutionArray[i]+k1*(2/3), timeArray[i]+deltat*(2./3.), **optional)[0]
        #Next integration step:
        SolutionArray[i + 1] = SolutionArray[i] + (k1/4. + (3./4.)*k2)
        
    return SolutionArray
