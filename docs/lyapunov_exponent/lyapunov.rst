Lyapunov exponent computation
=============================

The computation of the largest lyapunov exponent is based on a long time integration of two nearby trajectories, in which the separation in calculated at fixed time intervals via a recursive rescaling of the mutual distance between the two trajectories.

.. figure:: ../../img/lyapunov_rescaling.pdf
   :align: center
   :width: 65%

The following formalism applies.

The leading trajectory starts from an initial value :math:`\mathbf{x}_1`. Fixed a sampling time :math:`\Delta t`, the :math:`(i+1)-th` point of the reference trajectory is:

.. math:: 
   :nowrap:
 
   \begin{equation}
   \mathbf{x}_{i+1} = f\big(\mathbf{x}_{i}\big) \Big \rvert_{t_0}^{t_{0} + \Delta t} 
   \end{equation}


The neighbouring trajectory is initially perturbed with a perturbation :math:`\delta \mathbf{x}`. At every sampling point this trajectory is riscaled with the module of the initial perturbation at time :math:`t_0`.

The exact separation of the trajectories is:

.. math:: 
   :nowrap:
 
   \begin{equation}
   \mathbf{d}_{i+1} = f\big(\mathbf{z}_{i}\big) \Big \rvert_{t_0}^{t_{0} + (i+1)\Delta t}- \mathbf{x}_{i+1}
   \end{equation}

the vector :math:`\mathbf{d_{i+1}|}` is re-scaled with the module of the initial perturbation :math:`|\delta \mathbf{x}|`.

The next integration point for the perturbed trajectory is therefore:

.. math:: 
   :nowrap:
 
   \begin{equation}
   \mathbf{z}_{i+1} = \mathbf{x}_{i+1} + \mathbf{d}_{i+1} \frac{|\delta \mathbf{x}|}{|\mathbf{d}_{i+1}|}
   \end{equation}

By recursively finding the :math:`\mathbf{z}_{i}` points via rescaling the distance, and then measuring the deviation of the two trajectories, the estimation of the largest Lyapunov exponent is calculated then from the series:

.. math:: 
   :nowrap:
 
   \begin{equation}
   \lambda_{t} = \frac{1}{N \Delta t} \sum_{i=1}^{N} \ln \Big(\frac{\delta \mathbf{x}}{|\mathbf{d}_{i}|}\Big)
   \end{equation}