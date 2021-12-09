Rossler's system
================

Rossler's system is described by the system of differential equations:

.. math::
   \dot{x}&= -y - z\\
   \dot{y} &= x + ay \\
   \dot{z} &= b + z(x-c) \\

The stability matrix of this system is:

.. math::
   \begin{equation}
   \mathbb{A}(\mathbf{x}(t), t) =
   \begin{pmatrix}
   0 & -1 & -1\\
   1 & a & 0\\
   z & 0 & x-c
   \end{pmatrix}
   \end{equation}

1. Generate the input file containing the ODE system and the hard code the stability matrix, inside ``multiflap/odes/rossler.py``:

.. code-block:: python

   
        import numpy as np


        """
        Case adopted from:

            Optimized shooting method for finding periodic orbits of nonlinear dynamical systems
            Dednam, W and Botha, Andre E
            https://arxiv.org/abs/1405.5347

        """
        class Rossler:
            def __init__(self, a=None, b=None, c=None, d=None, e=None, q=None, p=None):

                self.a = a
                self.b = b
                self.c = c
                self.d = d
                self.e = e
                self.q = q
                self.p = p
                self.dimension = 3
            def dynamics(self, x0, t):

                """ODE system
                This function will be passed to the numerical integrator

                Inputs:
                    x0: initial values
                    t: time

                Outputs:
                    x_dot: velocity vector
                """
                x, y, z = x0
                dxdt = - y - z
                dydt = x + self.a*y
                dzdt = self.b + z*(x - self.c)
                vel_array = np.array([dxdt, dydt, dzdt], float)
                return vel_array



            def get_stability_matrix(self, x0, t):

                """
                Stability matrix of the ODE system

                Inputs:
                    x0: initial condition 
                Outputs:
                    A: Stability matrix evaluated at x0. (dxd) dimension
                    A[i, j] = dv[i]/dx[j]
                """
                x, y, z = x0
                A_matrix = np.array([[0, -1, -1],
                                    [1, self.a, 0],
                                    [z, 0, x-self.c]], float)

                return A_matrix

2. Generate the main file to run in the directory ``multiflap/main_rossler.py``:

Import the class generated in the input file ``Rossler`` and the modules to run and solve the multiple-shooting

.. code-block:: python

   from  odes.rossler import Rossler
   from ms_package.rk_integrator import rk4
   from ms_package.multiple_shooting_period import MultipleShooting
   from ms_package.lma_solver_period import Solver

set the initial guess:

.. code-block:: python
   
   x = [10., 10., 3.6]

Generate the object containing the Rossler's equations:

.. code-block:: python

   mymodel = Rossler(a=0.2, b=0.2, c=5.7)

Passe the object to the multiple-shooting class, and solve it

.. code-block:: python

   ms_obj =  MultipleShooting(x, M=2, period_guess= 5., t_steps=50000, model=mymodel)
   mysol = Solver(ms_obj = ms_obj).lma(5.)

.. toggle-header::
    :header: ``rossler_main.py`` **Show full main**

            .. code-block:: python

                import numpy as np
                from  odes.rossler import Rossler
                from ms_package.rk_integrator import rk4
                from ms_package.multiple_shooting_period import MultipleShooting
                from scipy.integrate import odeint
                import matplotlib.pyplot as plt
                from ms_package.lma_solver_period import Solver

                x = [10., 10., 3.6]

                time_array = np.linspace(0, 180, 90000)
                mymodel = Rossler(a=0.2, b=0.2, c=5.7)

                ms_obj =  MultipleShooting(x, M=2, period_guess= 5., t_steps=50000, model=mymodel)

                mysol = Solver(ms_obj = ms_obj).lma(5., '/Users/gducci/UCL/PROJECT/Simulations/class_test')

                jac = mysol[4]

                eigenvalues, eigenvectors = np.linalg.eig(jac)


                sol_array = mysol[3].space
                sol_time = mysol[3].time
                period = sol_time[-1]

                plt.plot( sol_time, sol_array[:,0], label = "D1")
                plt.plot( sol_time, sol_array[:,1], label = "D2")
                plt.plot( sol_time, sol_array[:,2], label = "R")
                plt.legend()
                plt.show()

The relocation of the points is shown in the following animation.

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/xuCNVc1aVBE" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

