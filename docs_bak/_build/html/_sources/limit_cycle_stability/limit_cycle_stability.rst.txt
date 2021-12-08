Stability of cycles
===================

Given the generic system of ordinary differential equations of the form:

.. math::    
   :label: odes
   
   \dot{\mathbf{x}} = \mathbf{v}(\mathbf{x}, t)

we want to find a particular solution such that:

.. math:: 
   :nowrap:
 
   \begin{equation}
   \mathbf{x}^{*}(t) = \mathbf{x}^{*}(t + T)
   \end{equation}

where :math:`T > 0` is the period of the cycle.

.. figure:: ../../img/limit_cycle_stability.png
   :align: center
   :alt: Limit cycle stability.
   :width: 65%

Depending on their asymptotic behavior limit cycles can be either **stable** or **unstable**. In the first scenario they behave as attractors, while in the second one they behave as repellers.

Let  :math:`\textbf{x}^*(t) + \delta\textbf{x}(t)` be the neighbor of a limit cycle state at time :math:`t`, where :math:`\delta\textbf{x}(t)` thus captures the initial perturbation with respect to the limit cycle. After a time equal to the period :math:`T`, this point of the state space is transported by the flow as following:

.. math::
   :label: flow_map

   x_{i}^{*}(t) + \delta{x}_{i}(t+T) = f_{i}\big(\textbf{x}^*(t) + \delta\textbf{x}(t)\big) \Big \rvert_{t}^{t+T}

By expanding the right hand side of Equation :eq:`flow_map` to the first order, and considering that :math:`x_{i}^{*}(t) = f_{i}\big(\textbf{x}^*(t)\big) \big \rvert_{t}^{t+T}` (limit cycle condition), we obtain:

.. math::
   :label: flow_map_2
  
   \delta{x}_{i}(t+T) = \sum_{j = 1}^{N} \frac{\partial f_{i}\big(\textbf{x}^*(t)\big) \big \rvert_{t}^{t+T}}{\partial x_{j}}\cdot \delta x_j(t)

Equation :eq:`flow_map_2` describes the dynamic of the perturbed neighboring state.

.. Note:: 
   The quantity

   .. math:: \mathbb{J}(\mathbf{x}^*)\big \rvert_{t}^{t+T}= \sum_{j = 1}^{N} \frac{\partial f_{i}\big(\textbf{x}^*(t)\big) \big \rvert_{t}^{t+T}}{\partial x_{j}}

   is called the Jacobian Matrix or Floquet Matrix, and describes how the neighbor is deformed by the flow, after a period :math:`T`.

The stability of a periodic orbit is governed by the eigenvalues of the Jacobian matrix (Floquet multipliers). For an orbit to be stable, all the eigenvalues (except the trivial one) must be less the one in absolute value.


The Jacobian (or Monodromy) matrix, is the solution of the differential equation:

.. math::
   :label: jacode

   \frac{d\mathbb{J}}{dt}(\mathbf{x}_0) \Big \rvert_{t_0}^{t}&= \mathbb{A}(\mathbf{x}, t) \mathbb{J}(\mathbf{x}_0) \Big \rvert_{t_0}^{t}\\
   \mathbb{J} (\mathbf{x}_0)\Big \rvert_{t_0}^{t_0} &= \mathbb{I}


.. Note:: 
   We define **stability matrix** the quantity  

   .. math:: \mathbb{A}(\mathbf{x}, t) = \frac{\partial}{\partial x_j}v_{i(\mathbf{x}, t)}
   
   where :math:`v` is the right-hand-side of Equation :eq:`odes`

Two methods have been implemented in the code to calculate the Jacobian matrix: 

* Analytical, solving Equation :eq:`jacode`;

* Numerical, by finite differentiation.

A **multiple-shooting** algorithm has been implemented in order to detect limit cycles and assess their stability.
