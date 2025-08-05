# xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
#
# Copyright (C) 2024 Dmitrii Kulik, S.Dmytriyeva
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

from xgems import *
from numpy import *

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

engine = ChemicalEngineDicts("resources/CalciteBC/CalciteBC-dat.lst")

engine.T = 298.15
engine.P = 100000.0
bulk_composition = {'C': 1e-08, 'Ca': 1e-08, 'Cl': 0.002, 'H': 111.016746657646,
                    'Mg': 0.001, 'O': 55.5083933588231, 'Sn': 130.841288437146, 'Zz': 0.0}
engine.set_bulk_composition(bulk_composition)
engine.supress_multiple_species(['ClO4-', 'Cl-'])

engine.equilibrate();

print("pH", engine.pH)
print("pE", engine.pE)
print("ionic_strength", engine.ionic_strength)
print("system_volume", engine.system_volume)
print("system_mass", engine.system_mass)

phase_names = engine.phase_names
out1 = engine.phase_sat_indices
out2 = engine.phases_molar_volume
out3 = engine.phases_moles
out4 = engine.phases_mass
out5 = engine.phases_volume
out6 = engine.phases_volume_frac
out7 = engine.solids_mass_frac

#phases
print("\n           phase; phase_sat_indices; phases_molar_volume;  phases_moles;     phases_mass;   phases_volume; phases_volume_frac; solids_mass_frac;")
for name in phase_names:
    print("{:16s}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g};".format(
           name, out1[name], out2[name], out3[name], out4[name], out5[name], out6[name], out7[name]))

#species
species_names = engine.species_names
out1 = engine.species_charges
out2 = engine.species_molar_mass
out3 = engine.species_molar_volumes
out4 = engine.species_moles
out5 = engine.species_ln_activities
out6 = engine.species_ln_activity_coefficients

print("\n      specie;  charges;      molar_mass;     molar_volumes;            moles;    ln_activities; ln_activity_coefficients;")
for name in species_names:
    print("{:12s}; {:8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g};".format(
           name, out1[name], out2[name], out3[name], out4[name], out5[name], out6[name]))
