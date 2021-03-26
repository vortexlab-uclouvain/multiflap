[![Static](https://img.shields.io/badge/docs-latest-blue)][stable docs url]

[stable docs url]: https://multiflap.readthedocs.io/en/latest/
# multiflap


`multiflap` is a Python toolbox for finding periodic orbit of nonlinear systems
of ODEs. The code relies on the multiple-shooting algorithm, and simoultaneously allows to study the stability of the system via the Floquet multipliers.

It handles also the computation of the leading Lyapunov exponent to characterize chaotic behaviours. This computation relies on long time numerical calculation of two nearby trajectories, rescaling the separation of the perturbed one in defined sampling points.

This toolbox is optimised to study the limit cycle of flapping wings dynamics, but the modular architecture allows the user to use the kernel with any set of ODEs.

In order to run the code for a particular dynamical system, two files have to be created:
* **Input file**: The input file contains the system of differential equations, and the related staility matrix.
* **Main file**: The main file is the python module that will e executed. More detailes on the format of the main file will be provided later. **NOTE** despite of example cases, the main file has to be inside `multiflap` directory.

For guided examples and tutorials of how to use multiflap, and how to create your case-study, please refer to the dedicated documentation available [here](<https://multiflap.readthedocs.io/en/latest/tutorial/index.html>).

# Contacts
Gianmarco Ducci, PhD student @UCLouvain, Belgium.

<gianmarco.ducci@uclouvain.be>
## Installation and getting started

`multiflap` runs on Python 3 and Python 2.  

1.   Get a local copy of `multiflap`:

```
git clone https://github.com/vortexlab-uclouvain/multiflap.git && cd multiflap
```
2. Install the packages:
```bash
python setup.py install
```
3. To run one of the examples, from the multiflap root directory:

```
python3 examples/rossler_system.py
```
## Input files

The input file is the file containing the set of ODEs, and saved within the directory
```
multiflap/odes
```

Examples of input file can be found [here](multiflap/odes).
The name of the class is arbitrary, while the names of the methods (`dynamics` and `get_stability_matrix`) cannot be changed because this is the interface that communicates with the multiple-shooting kernel. If the stability matrix of the system is hard to hand-code, the user can run the code by setting the numerical computation from the aguments of `MultipleShooting` or `MultipleShootingPeriod` classes, as shown in the [example](examples/lorentz_system.py).

## Main file

Two different classes should be called, depending if the period of the periodic orbit is known or not a priori.

#### Unknown period

For periodic orbits with unknown periods, the code relies on the modules `multiple_shooting_period.py` and `lma_solver_period.py`.

1. Import the model module `your_model.py`, and the other modules needed to run and solve the multiple shooting problem. Then create the model object. The model is created by calling the Class previously defined within [odes](multiflap/odes) directory.

```python
from odes.your_model import YourOdeClass
from ms_package.multiple_shooting_period import MultipleShootingPeriod
from ms_package.lma_solver_period import SolverPeriod
your_model = YourOdeClass() 	# takes eventual arguments (parameters) if needed
```
2. Call the multiple-shooting scheme. It take as arguments:

* **x**: state space points guess
* **M**: int, number of multiple-shooting points
* **period_guess**: first guess of the periodic orbit period
* **t_steps**: time steps from two consecutive points
* **model**: ODE class
* **option_jacobian**: 'numerical' or 'analytical' (if not specified, analytical is the default)

```python
ms_obj = MultipleShootingPeriod(x, M=2, period_guess= 23.,
		t_steps=50000, model=your_model, option_jacobian = 'analytical')
```
3. Call the multiple-shooting solver.
```python
mysolution = SolverPeriod(ms_obj = ms_obj).lma()
```

#### Known period 


For periodic orbits with known periods, the code relies on the modules `multiple_shooting.py` and `lma_solver.py`.

1. Build the model object. The model is created by calling the Class previously defined within [odes](multiflap/odes) directory.

```python
from odes.your_model import YourOdeClass
from ms_package.multipleshooting import MultipleShooting
from ms_package.lma_solver import Solver
your_model = YourOdeClass()
```
2. Call the multiple-shooting scheme. It take as arguments:

* **x**: state space points guess
* **M**: int, number of multiple-shooting points
* **period**: periodic orbit period
* **t_steps**: time steps from two consecutive points
* **model**: ODE class
* **option_jacobian**: 'numerical' or 'analytical' (if not specified, analytical is the default)

```python
ms_obj = MultipleShooting(x, M=2, period= 0.25,
		t_steps=50, model=your_model, option_jacobian = 'analytical')
```
3. Call the multiple-shooting solver.
```python
mysolution = Solver(ms_obj = ms_obj).lma()
```
## Generate a phase portrait animation

A script that allows to produce phase portrait animation using `Matplotlib` is provided within the [utilities](multiflap/utils) directory. In the default template the system of ODE is given, however is highly recommended to use this script as mere post-processing and load data from simulations as arrays.
Outputs of the scripts, as examples, are within [anim](anim) directory.


## Literature

#### Multiple shooting scheme 

Improved numerical Floquet multipliers \[[pdf](pdfs/lust2001.pdf)\]\
Lust, Kurt\
International Journal of Bifurcation and Chaos, vol. 11, pp. 2389--2410, 2001

ChaosBook, Chaos: Classical and Quantum \
P. Cvitanovi\'c, R. Artuso, R. Mainieri, G. Tanner and G. Vattay \
Niels Bohr Institute, Copenhagen 2016

#### Example cases and validation

###### Rossler's system \[[odes](multiflap/odes/rossler.py)\]
Optimized shooting method for finding periodic orbits of nonlinear dynamical systems \[[pdf](pdfs/dednam2015.pdf)\]\
Dednam, W and Botha, Andre E \
Engineering with Computers, vol. 31, num. 4, pp. 749--762, 2015, Springer

Contruction of Poincaré Return Maps for Rössler Flow \[[pdf](pdfs/basu2007.pdf)\] \
A. Basu \
Tech. Report, School of Engineering, Georgia Institute of Technology, 2007 

###### Redox oscillator \[[odes](multiflap/odes/redox_oscillation.py)\]
A Robust Model for Circadian Redox Oscillations \[[pdf](pdfs/delolmo2019.pdf)\]\
del Olmo, M.; Kramer, A.; Herzel, H.
*Int. J. Mol. Sci.*, **2019**, 20, 2368

###### Forced Van der Pol oscillator \[[odes](multiflap/odes/forced_vdp.py)\]
Practical Bifurcation and Stability Analysis, pp. 347--349  
Seydel, R., 2009
###### Isothermal reaction \[[odes](multiflap/odes/isothermal_reaction.py)\]
Practical Bifurcation and Stability Analysis, pp. 325--327  
Seydel, R., 2009
