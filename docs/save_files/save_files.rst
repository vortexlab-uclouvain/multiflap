Save results of a simulation
=============================

To save the simulation results, the package comes with the class ``SaveData``.

In order to store your data, you first have to create the folder ``results`` on the top of the multiflap directory.

1. Creation of ``results`` directory

.. code-block:: bash

   $ mkdir results
 
Then inside the main file create the study case. This will automatically generate a subfolder within ``results/my_study_case`` in which the results will be stored and then save your data.

2. Inside the main file:

.. code-block:: python

       from utilities.save_data import SaveData

       sim_name = "my_study_case"
       SaveData(sim_name).make_folder() # creation of my_study_case subfolder within results
       # 
       # ... simulation
       #

       # save data

       SaveData(sim_name).save_data('my_results_name', my_array)

Example of save files in :ref:`Redox oscillation` (Show full main).
