Jacobian computation
====================

Two methods are implemented to calculate the Jacobian matrix and build the diagonal blocks of the multiple-shooting matrix :math:`\textbf{M}`.

Analytical computation
**********************

The analytical computation solves the Equation:


.. math::
   :label: jac_an

   \frac{d\mathbb{J}}{dt}(\mathbf{x}_0) \Big \rvert_{t_0}^{t}&= \mathbb{A}(\mathbf{x}, t) \mathbb{J}(\mathbf{x}_0) \Big \rvert_{t_0}^{t}\\
   \mathbb{J} (\mathbf{x}_0)\Big \rvert_{t_0}^{t_0} &= \mathbb{I}

This means solving a :math:`(n+n^{2})` system of differential equations.

The system is therefore accordingly re-shaped.

.. py:method::
Numerical computation
*********************
