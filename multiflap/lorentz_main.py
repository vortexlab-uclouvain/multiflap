'''
file lorentz_main.py
@author Gianmarco Ducci
@copyright Copyright © UCLouvain 2020

multiflap is a Python tool for finding periodic orbits and assess their stability via the Floquet multipliers.

Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>

List of the contributors to the development of multiflap, Description and complete License: see LICENSE and NOTICE files.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
import numpy as np
from  odes.lorentz import Lorentz
from ms_package.lyapunov_exponents import LyapunovExponents

x = [10., 10., 3.6]

# Creating the object containing the eqs
mymodel = Lorentz(a=10, b=28, c=8/3)

t_f = 100
n = 10000

# calling the LyapunovExponents class
lambda_t = LyapunovExponents(x, n, t_f, mymodel).get_lyapunov_exponent()

print(lambda_t)
