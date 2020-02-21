import numpy as np  # Import numpy
import FlappingForcesDev as FF
import scipy.integrate as ode

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

def RK2(velocityFunction, initialCondition, timeArray, **optional):
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

if __name__ == "__main__":
    #This block will be evaluated if this script is called as the main routine
    #and will be ignored if this file is imported from another script.
        #
        #This is a handy structure in Python which lets us test the functions
    
        #In order to test our integration routine, we are going to define Harmonic
        #Oscillator equations in a 2D state space:
        def birdEqn_py(t, ssp):
                """
                State space velocity function for the Equation of Motion of the longitudinal plane
        
                Inputs:
                ssp: State space vector
                ssp = (u, w, q, theta)
                t: Time
                
                Outputs:
                vel: Time derivative of ssp.
                """
                #Parameters:
                g = 9.81
                mass = 1.2
                #Read inputs:
                u, w, q, theta  = ssp  # Read state space points
                u = ssp[0]
                w = ssp[1]
                q = ssp[2]
                theta = ssp[3]
                # Longitudinal equation of motion:
                [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FF.FlappingForces(t, u, w, q, theta)
                dudt = -q*w - g*np.sin(theta) - Fz/mass 
                dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
                dqdt =  My/0.1
                dthetadt = q                # Collect Equations of motion in a single NumPy array:
        
                vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector            
                return vel


        def birdEqn(ssp, t):
                """
                State space velocity function for the Equation of Motion of the longitudinal plane
        
                Inputs:
                ssp: State space vector
                ssp = (u, w, q, theta)
                t: Time
                
                Outputs:
                vel: Time derivative of ssp.
                """
                #Parameters:
                g = 9.81
                mass = 1.2
                #Read inputs:
                u, w, q, theta  = ssp  # Read state space points
                u = ssp[0]
                w = ssp[1]
                q = ssp[2]
                theta = ssp[3]
                # Longitudinal equation of motion:
                [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FF.FlappingForces(t, u, w, q, theta)
                dudt = -q*w - g*np.sin(theta) - Fz/mass
                dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
                dqdt =  My/0.1                               #(0.01*Fy - 0.1*(Fz))/0.1 # <- Reference
                dthetadt = q
                # Collect Equations of motion in a single NumPy array:
        
                vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector            
                return vel

        tInitial = 0
        tFinal = 0.25
        Nt = 60  # Number of points time points in the interval tInitial, tFinal
        tArray = np.linspace(tInitial, tFinal, Nt)
    
        #Initial condition for the Harmonic oscillator:
        ssp0 = np.array([1.0, 0], float)
        ssp0_eqn = np.array([17.89960883175204, -1.4086561216686262, 0.022401202083810585, -0.1678545149865401], float)

        #Compute the solution using Runge-Kutta routine:
        sspSolution = RK4(birdEqn, ssp0_eqn, tArray)
        print("RK done")
        
        sol = ode.solve_ivp(birdEqn_py, [tInitial, tFinal], ssp0_eqn,'BDF', tArray)
        print("Python routine done")
        
        #from scipy.integrate import odeint
    #        sspSolution = odeint(velocity, ssp0, tArray)
        xSolution = sspSolution[:, 0]
        vSolution = sspSolution[:, 1]
        import matplotlib.pyplot as plt
        
        plt.plot(sol.t, sol.y[0])
        plt.scatter(tArray, xSolution, color = 'red')
        plt.grid(True)
        plt.show()
    
        
