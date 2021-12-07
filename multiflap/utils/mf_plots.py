import numpy as np
import matplotlib.pyplot as plt

class Plot:

    def __init__(self, color='k', grid=True):
        self.color = color
        self.grid = grid
    
    def limit_cycle_2D(self, x, y):
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(x, y, color=self.color)
        plt.grid(self.grid)
        plt.show()
        return None

    def limit_cycle_3D(self, *args):

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