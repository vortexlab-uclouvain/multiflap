# multiflap

`multiflap` is a Python toolbox for finding periodic orbit of nonlinear systems
of ODEs. The code relies on the multiple-shooting algorithm, and simoultaneously allows to study the stability of the system via the Floquet multipliers.

This toolbox is optimised to study the limit cycle of flapping wings dynamics, but the modular architecture allows the user to use the kernel with any set of ODEs, changing an input file only.

## Installation and getting started

`multiflap` runs on Python 3 and Python 2.  

1.   Get a local copy of `multiflap`:

```
git clone git@git.immc.ucl.ac.be:gducci/multishooting_flapping.git 
```
2. Add multiflap to the PYTHONPATH environment variable

```
bash set_evironment_variable
```
3. To run one of the examples, from the multiflap root directory:

```
python3 examples/rossler_system.py
```
## Input files

The input file is the file containing the set of ODEs. The name of the class is arbitrary, while the names of the methods cannot be changed because this is the interface that communicates with the multiple-shooting kernel.

```
class 'MyClass'
	def __init__(self, a=None, b=None, c=None, d=None, e=None, q=None, p=None):

		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.e = e
		self.q = q
		self.p = p
	
	def 'dynamics'(self, x0, t):

		# type here ODEs
		return vel_array
	
	def 'get_stability_matrix'(self, x0, t):

		# type here the stability matrix

		return A_matrix

```