import numpy as np
import matplotlib.pyplot as plt

class Plot():

    """
    Plotting module for quick visualization of simulation results. This class
    is not intended to have customized and high-quality images for publication,
     rather a quick outcome of simulations.

    Example:
        >>> import multiflap as mf
        >>> plot = mf.Plot(color='red') # create the plot object
        >>> plot.limit_cycle_2D(sol_array[:,0], sol_array[:,1])
        
    """

    def __init__(self, color='k', grid=True):
        self.color = color
        self.grid = grid

    def limit_cycle_2D(self, x, y):
        """Plot of 2D limit cycles in the phase space.

        Parameters
        ----------
        x : 1D-array
            Array of the first variable to be plotted
        y : 1D-array
            Array of the second variable to be plotted

        Returns
        -------
        None : void function
               It plots the phase space of the 2 variables


        """
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(x, y, color=self.color)
        plt.grid(self.grid)
        plt.show()
        return None

    def limit_cycle_3D(self, *args):
        """Plot of 3D limit cycles in the phase space.

        Parameters
        ----------
        args : [nxm] array
               Array solution, output of multiflap

        Returns
        -------
        None : void function
               It plots the phase space of the first 3 variables of args


        """

        if len(args) == 2:
            raise ValueError('Please use the method "limit_cycle_2D"')
        x = args[0]
        y = args[1]
        z = args[2]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$') 
        ax.set_zlabel('$z$') 
        ax.plot(x,
                y,
                z)
        plt.show()

        return None

    def plot_multipliers(self, eigenvalues):
        """Plot the Floquet multipliers and the unitary circle

        Parameters
        ----------
        eigenvalues : array
                      eigenvalues array, output of ``multiflap``

        Returns
        -------
        None : void function
               It plots the Floquet multipliers with the unitary circle


        """
        # Define the unitary circle
        circle = np.linspace(0,2*np.pi,101)
        fig = plt.figure()
        ax = fig.gca()
        ax.set_aspect('equal')
        ax.set_xlabel('$Re$', fontsize=20)
        ax.set_ylabel('$Im$', fontsize=20)
        ax.scatter(eigenvalues.real, eigenvalues.imag, color=self.color)
        ax.plot(np.cos(circle), np.sin(circle), color='cornflowerblue')
        plt.grid(self.grid)
        plt.show()
        return None

    def plot_time_series(self,
                         variables,
                         time,
                         ylabel='variable',
                         label=False):

        # counting the variables to be plotted 
        var_number = len(variables[0,:]) 

        fig= plt.figure()
        ax = fig.gca()
        ax.set_xlabel('time [s]')
        ax.set_ylabel(ylabel)
        for i in range(var_number):
            if label!= False:
                ax.plot(time, variables[:,i], label=label[i])
            else:
                ax.plot(time, variables[:,i])
        plt.legend()
        plt.show() 
        return None
