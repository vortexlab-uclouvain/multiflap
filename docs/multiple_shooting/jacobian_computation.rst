Jacobian computation
====================

Two methods are implemented to calculate the Jacobian matrix and build the diagonal blocks of the multiple-shooting matrix :math:`\textbf{M}`.

Analytical computation
**********************

The analytical computation solves the Equation:

.. math::
   :label: jac_an

   \begin{pmatrix} \dot{\mathbf{x}}\\ \dot{\mathbb{J}} \end{pmatrix} = \begin{pmatrix} \mathbf{v}(\mathbf{x}, t)\\ \mathbb{A}(\mathbf{x}, t) \ \mathbb{J} \end{pmatrix}

with the initial conditions:

.. math::
   :label: jac_an_initial_conditions

   \begin{pmatrix} {\mathbf{x}(t_0)}\\ {\mathbb{J}^{0}} \end{pmatrix} = \begin{pmatrix} \mathbf{x}_0 \\ \mathbb{I} \end{pmatrix}

This means solving a :math:`(n+n^{2})` system of differential equations.

First of all the ODE system :eq:`jac_an` is reshaped in a :math:`(n+n^{2})` equations. This is done by the method ``jacobian_ode``.

.. toggle-header::
    :header: ``jacobian_ode`` **Show code**

            .. code-block:: python

                    def jacobian_ode(self, x0_jacobian, t):
                        """
                        Set up the additional ODE system (d + d^2) to evaluate analytically the
                        Jacobian matrix.
                        This function is used to unpack the Jacobian elements to solve

                        \dot{J} = AJ

                        It reshapes Equation (7.18) of Seydel's book "Practical Bifurcation and
                        Stability Analysis"
                        in order get the components of the Jacobian via numerical integration
                        Inputs:
                            x0_jacobian: (d+d^2) initial values
                                         state space itself and the tangent space
                            t: time.

                        Outputs:
                            jac_elements_ODE = (d+d^2) dimensional velocity vector

                        """



                        x0 = x0_jacobian[0:self.dim]
                        J = x0_jacobian[self.dim:].reshape((self.dim, self.dim))

                        jac_elements_ODE = np.zeros(np.size(x0_jacobian))

                        jac_elements_ODE[0:self.dim] = self.model.dynamics(x0, t)


                        #Last dxd elements of the velJ are determined by the action of
                        #stability matrix on the current value of the Jacobian:

                        velocity_vector = np.dot(self.model.get_stability_matrix(x0, t), J)

                        # shape a back a (dxd) array in a d^2 matrix

                        jac_elements_ODE[self.dim:] = np.reshape(velocity_vector,
                                                                       self.dim**2)

                        return jac_elements_ODE

The function ``get_jacobian_analytical`` is then implemented in order to solve the reshaped ODE system :eq:`jac_an`

.. toggle-header::
    :header: ``get_jacobian_analytical()`` **Show code**

            .. code-block:: python

                    def get_jacobian_analytical(self, x0, initial_time, integration_time):

                        """
                        Return the Jacobian (or Monodromy Matrix) of the flow, starting from x0
                        and integrated for a time "integration_time".

                        It solves numerically the (d + d^2) ODE system.

                        Reference: Equation (7.18) Seydel's book "Practical Bifurcation and
                        Stability Analysis".

                        Inputs:
                            x0 : Initial point of the phase space. len(x0) = dimension
                            initial_time: initial time needed as the system is non-autonomous
                            integration_time: integration time
                        Outputs:
                            J: Jacobian (Monodromy Matrix) of flow from t -> t+integration_time
                        """

                        # Initial conditions (ic) for Jacobian matrix (see 7.18 Seydel)

                        jac_ic = np.identity(self.dim)

                        jacODE_ic = np.zeros(self.dim + self.dim ** 2)

                        jacODE_ic[0:self.dim] = x0

                        jacODE_ic[self.dim:] = np.reshape(jac_ic, self.dim**2)


                        t_Final = initial_time + integration_time
                        Nt = 50000 #50000  # interval discretization for computing the integration

                        tArray = np.linspace(initial_time, t_Final, Nt)

                        start_jac = time.time()

                    #   jac_elements_solution = ode.solve_ivp(jacobian_ode,[t_initial, t_Final],
                                                #jacODE_ic, 'RK45')
                        if self.integrator=='rk':
                            rk_jac_elements_solution = rk2(self.jacobian_ode, jacODE_ic, tArray)

                            end_jac = time.time()
                            print("Jacobian time ", (end_jac-start_jac))

                            jac_elements_solution = rk_jac_elements_solution.x
                        #    jac_elements_solution = jac_elements_solution.y.T
                            # Pack back the jacobian in matrix:
                            J = jac_elements_solution[-1, self.dim:].reshape((self.dim,
                                                                                self.dim))
                        if self.integrator=='odeint':
                            jac_elements_solution = odeint(self.jacobian_ode, jacODE_ic, tArray)
                            J = jac_elements_solution[-1, self.dim:].reshape((self.dim,
                                                                              self.dim))
                        return J
Numerical computation
*********************

Alternatively, Jacobian matrix can be computed via numerical differentiation of the perturbed
trajectory along each state variable, i.e.:

.. math::
  :label: jac_numerical 

  \mathbb{J}_{i,j} (\mathbf{x}_0) \Big \rvert_{t}^{t+T} = \frac{f_{i}(\mathbf{x_0} + \varepsilon\mathbf{\hat{e}}_{j}) \Big \rvert_{t}^{t+T} -f_{i}(\mathbf{x_0}) \Big \rvert_{t}^{t+T} }{\varepsilon} 

.. tip::

   By computing the Jacobian numerically, there is no need to hard code the **stability matrix**. 
   This feature is particularly useful anytime that the encoding of the stability matrix is hard to provide.

The function that returns the numerical Jacobian is ``get_jacobian_numerical``.

.. toggle-header::
    :header: ``get_jacobian_numerical`` **Show code**

            .. code-block:: python

                    def get_jacobian_numerical(self, x0, initial_time, integration_time):

                        """
                        Finite difference evaluation of Jacobian
                        Flow is perturbed in all directions of the phase space, and the generic
                        component of the Jacobian is calculated by finite difference as follow:

                        dF[i]/dx[j] = (F^t(x_perturbed)[i] - F^t(x)[i])/perturbation

                        Jacobian = dF[i]/dx[j] (dim x dim) matrix

                        epsilon = value of the perturbation
                        """
                        time_steps = 50000
                        # ---------------------------------------------------------------------
                        #  Initialization of the Jacobian Matrix
                        # ---------------------------------------------------------------------

                        jacobian = np.zeros((self.dim,self.dim))

                        # ---------------------------------------------------------------------
                        # Set the numerical perturbation over the direction of the flow
                        # ---------------------------------------------------------------------

                        epsilon = 1e-3


                        # ---------------------------------------------------------------------
                        # Finite difference scheme for jacobian evaluation
                        # ---------------------------------------------------------------------

                        print("... Running jacobian Function")
                        for j in range (self.dim):
                            perturbation = np.zeros(self.dim)

                            perturbation[j] = perturbation[j] + x0[j]*epsilon

                            x0_pert = x0 + perturbation

                            [vel, _] = self.get_mappedpoint(x0,
                                                            initial_time,
                                                            integration_time)

                            [vel_pert, _] =  self.get_mappedpoint(x0_pert,
                                                                  initial_time,
                                                                  integration_time)


                            for i in range (self.dim):

                                jacobian[i,j] = (vel_pert[i] - vel[i])/perturbation[j]


                        print("... jacobian Calculated")
                        return jacobian
                               
                   
