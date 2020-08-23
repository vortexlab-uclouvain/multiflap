Multiple-shooting algorithm
===========================

What follows is the description of the **multiple-shooting** algorithm. Two functions are implemented, depending whther the orbit period is known a priori (i.e. for non-autonomous systems) or not (i.e. most of the autonomous systems). 

Both methods are adaptations of the algorithms described by Lust [cite here]. 

Two methods to compute the Jacobian (or Monodromy matrix) are also described.

.. toctree::
   :maxdepth: 2

   numerical_scheme.rst
   jacobian_computation.rst
   floquet_multipliers.rst
