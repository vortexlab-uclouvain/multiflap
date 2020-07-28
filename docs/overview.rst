Overview
========

`multiflap` is a Python toolbox for finding periodic orbit of nonlinear systems of ODEs. The code relies on the multiple-shooting algorithm, and simoultaneously allows to study the stability of the system via the Floquet multipliers.

This toolbox is optimised to study the limit cycle of flapping wings dynamics, but the modular architecture allows the user to use the kernel with any set of ODEs.

In order to run the code for a particular dynamical system, two files have to be created:
* **Input file**: The input file contains the system of differential equations, and the related staility matrix.
* **Main file**: The main file is the python module that will e executed. More detailes on the format of the main file will be provided later. **NOTE** despite of example cases, the main file has to be inside `multiflap` directory.
## Installation and getting started

`multiflap` runs on Python 3 and Python 2.  

1.   Get a local copy of `multiflap`:

```
git clone https://github.com/vortexlab-uclouvain/multiflap.git
```
2. Install the packages:
```bash
python setup.py install
```
3. To run one of the examples, from the multiflap root directory:

```
python3 examples/rossler_system.py
```
