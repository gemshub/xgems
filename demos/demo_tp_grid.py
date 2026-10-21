# xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
#
# Copyright (C) 2018 Allan Leal, Dmitrii Kulik
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import xgems as xg

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

xg_system_filename = 'resources/gems_files/M1-6-dat.lst'

xgEngine = xg.ChemicalEngine(xg_system_filename)

# preliminary check of thermodynamic consistency
code = xgEngine.reequilibrate()
print('GEMS return code:', code)

# element indexes
element_indexes = {
    'Pb': xgEngine.indexElement('Pb'),
    'I':  xgEngine.indexElement('I'),
    'O':  xgEngine.indexElement('O'),
    'He': xgEngine.indexElement('He'),
}

# species indexes
species_indexes = {
    'Pb(s)':   xgEngine.indexSpecies('Pb'),
    'I2(s)':   xgEngine.indexSpecies('I2'),
    'I(g)':    xgEngine.indexSpecies('I(g)'),
    'I2(g)':   xgEngine.indexSpecies('I2(g)'),
    'Pb(g)':   xgEngine.indexSpecies('Pb(g)'),
    'Pb2(g)':  xgEngine.indexSpecies('Pb2(g)'),
    'PbI(g)':  xgEngine.indexSpecies('PbI(g)'),
    'PbI2(g)': xgEngine.indexSpecies('PbI2(g)'),
    'PbI2(s)': xgEngine.indexSpecies('PbI2(s)'),
}

print('Element indexes:', element_indexes)
print('Species indexes:', species_indexes)

b0 = xgEngine.elementAmounts().copy()

b0[element_indexes['Pb']] = 1.61e-5
b0[element_indexes['I']]  = 2 * b0[element_indexes['Pb']]
b0[element_indexes['He']] = 0.75
b0[element_indexes['O']]  = 2e-10

P = 1e5  # Pa

He_amount = b0[element_indexes['He']]
print('He amount: ', He_amount, ' mol, P: ',P,' Pa')

temperatures = np.arange(25, 1000, 1)  # °C

results = {'T (°C)': [], 'status': []}
for name in species_indexes:
    results[name] = []

for T in temperatures:
    code = xgEngine.equilibrate(T + 273.15, P, b0)

    results['T (°C)'].append(T)
    results['status'].append(code)

    for name, idx in species_indexes.items():
        results[name].append(xgEngine.speciesAmount(idx))

    if code not in (2, 6):
        print('Warning: equilibrium not fully converged at {T:.1f} °C (code {code})')

print('Calculation finished.')
print(results)




