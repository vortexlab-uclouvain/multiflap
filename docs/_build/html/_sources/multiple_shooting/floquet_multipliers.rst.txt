Floquet multipliers
===================

Once the convergence is reached the Floquet multipliers can be automatically exctracted from the output of LMA solver.

The global Jacobian over one period is calculated considering the semi-group property:

.. math::
   :label: jacobian_period

   \mathbb{J} (\mathbf{x}_0) \Big \rvert_{0}^{T} =    \mathbb{J} (\mathbf{x}_{m-1}) \Big \rvert_{t_{m-1}}^{T} \cdots  \mathbb{J} (\mathbf{x}_1) \Big \rvert_{t_{1}}^{t_{1}+\tau}\cdot \mathbb{J} (\mathbf{x}_0) \Big \rvert_{0}^{\tau}
   
Eigenvalues (Floquet multipliers) and eigenvectors, are then calculated with ``python`` builtin functions from ``linalg``.
