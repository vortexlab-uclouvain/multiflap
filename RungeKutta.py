import numpy as np  # Import numpy
import flapping_forces as FF
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


if __name__ == "__main__":

    #This block will be evaluated if this script is called as the main routine
    #and will be ignored if this file is imported from another script.
        #
        #This is a handy structure in Python which lets us test the functions
    
        #In order to test our integration routine, we are going to define Harmonic
        #Oscillator equations in a 2D state space:
        lift_py = []
        moment_py = []
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
                [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FF.flapping_forces(t, u, w, q, theta)
                lift_py.append(Fy)
                moment_py.append(My)
                dudt = -q*w - g*np.sin(theta) - Fz/mass 
                dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
                dqdt =  My/0.1
                dthetadt = q                # Collect Equations of motion in a single NumPy array:
        
                vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector            
                return vel

        def birdEqn(ssp, t, **kinematics):
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
                [Fx, Fy, Fz, My, F_tail, M_wing, M_tail, M_drag, M_lift] = FF.FlappingForces(t, u, w, q, theta, **kinematics)
                lift.append(Fy)
                moment.append(My)
                dudt = -q*w - g*np.sin(theta) - Fz/mass
                dwdt = q*u + g*np.cos(theta) - Fy/mass - F_tail/mass
                dqdt =  My/0.1                               #(0.01*Fy - 0.1*(Fz))/0.1 # <- Reference
                dthetadt = q
                # Collect Equations of motion in a single NumPy array:
        
                vel = np.array([dudt, dwdt, dqdt, dthetadt], float)  # Velocity vector
                
                return vel
        lift = []
        moment = []

        tInitial = 0
        tFinal = .25
        Nt = 50  # Number of points time points in the interval tInitial, tFinal
        tArray = np.linspace(tInitial, tFinal, Nt)
    
        #Initial condition for the Harmonic oscillator:
        ssp0 = np.array([18., 0., 0., 0.], float)
        ssp0_eqn = np.array([18.5, -2.15, -0.11, -0.11], float)
        
        #Compute the solution using Runge-Kutta routine:
#        sspSolution = ode.odeint(birdEqn, ssp0, tArray)
        sspSolution_rk2 = ode.odeint(birdEqn, ssp0, tArray)
        sspSolution_rk = RK4(birdEqn, ssp0, tArray)

        print("RK done")
        
#        sol = ode.solve_ivp(birdEqn_py, [tInitial, tFinal], ssp0,'BDF', tArray)
#        print("Python routine done")
        lift = np.asanyarray(lift)
        moment = np.asanyarray(moment)
        lift_py=np.asanyarray(lift_py)
        moment_py = np.asanyarray(moment_py)
        lift_filtered = np.zeros(len(tArray))
        moment_filtered = np.zeros(len(tArray))
#        lift_filtered[0:-1] = lift[::2]
#        lift_filtered[-1] = lift[-1]
#        moment_filtered[0:-1] = moment[::2]
        moment_filtered[-1] = moment[-1]
#        xSolution = sspSolution[:, 0]
#        vSolution = sspSolution[:, 1]
        xSolution_RK4 = sspSolution_rk[:, 0]
        ySolution_RK4 = sspSolution_rk[:, 1]
        qSolution_RK4 = sspSolution_rk[:, 2]
        xthetaSolution_RK4 = sspSolution_rk[:, 3]
        
        xSolution_RK2 = sspSolution_rk2[:, 0]
        ySolution_RK2 = sspSolution_rk2[:, 1]
        qSolution_RK2 = sspSolution_rk2[:, 2]
        xthetaSolution_RK2 = sspSolution_rk2[:, 3]

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        Lift_LL = np.load('../utilities/forces_folder/lift_LL.npy')
        Moment_LL = np.load('../utilities/forces_folder/moment_LL.npy')
        LL_time = np.linspace(0., 0.25, len(Lift_LL))
#        plt.plot(sol.t, sol.y[0], label='solve_ivp')
        plt.plot(tArray, xSolution_RK4)
        plt.scatter(tArray, xSolution_RK2, color = 'red')

#        plt.scatter(tArray, lift_filtered, s=10, color='red', label='RK4')
#        plt.plot(LL_time, Lift_LL)
        plt.xlabel('t')
        plt.ylabel('$u(t)$')
        plt.grid(True)
        plt.show()
        plt.legend()
##        plt.savefig('/Users/gducci/UCL/MyWrittenDocs/Report/ValidationCoupling/figures/'+'flightdynamics.eps', format='eps')
#        fig15 = plt.figure()
#        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1.,1.])
#        
#        ax15 = fig15.add_subplot(gs[0, 0])
#        ax15.plot(tArray, lift_filtered,linewidth=2.0, color='blue', label='lift coupled')
#        ax15.plot(LL_time, Lift_LL  ,linewidth=2., color='red', label='lift uncoupled')
#        ax15.set(xlim=(0, 0.25))
#        ax15.set_ylabel('$L(t)$', fontsize=18)
#        plt.tick_params(labelsize=14)
#        plt.legend()
#        ax15.grid()
#        
#        ax15 = fig15.add_subplot(gs[1, 0])
#        ax15.plot(tArray, moment_filtered ,linewidth=2.0, color='blue', label='moment coupled')
#        ax15.plot(LL_time, Moment_LL  ,linewidth=2.,color = 'red', label='moment uncoupled')
#        ax15.set(xlim=(0, 0.25))
#        ax15.set_ylabel('$M_y(t)$', fontsize=18)
#        plt.tick_params(labelsize=14)
#        ax15.set_xlabel('$t$', fontsize=18)
#        plt.tick_params(labelsize=14)
#        plt.legend()
#        ax15.grid()
##        plt.savefig('/Users/gducci/UCL/MyWrittenDocs/Report/ValidationCoupling/figures/'+'coupling_forces.eps', format='eps')
#
#
#
#        
