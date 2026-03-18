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

from xgems import *
from numpy import *

#connect all loggers
update_loggers(True, "test_demo1.log", 3)

engine = ChemicalEngine("resources/solvus2/series1-dat.lst")

T = engine.temperature()
P = engine.pressure()
b = engine.elementAmounts()

engine.equilibrate(T, P, b)

# SetPT test
print("setPT(101325, 298.15)",  engine.setPT(101325, 298.15) )
print("setPT(101325, 673.15)",  engine.setPT(50000000, 673.15) )

# long names
ispecies1 = engine.indexSpecies("Anorthite_long_substance_name")
ispecies2 = engine.indexSpecies("Anorthite_long_substance_name", "Plagioclase")
ispecies3 = engine.indexSpecies("Anorthite")
print(ispecies1, " ",ispecies2, " ",ispecies3, " ")

print(engine.speciesName(ispecies1), " ", engine.speciesName(ispecies2), " ",engine.speciesName(ispecies3))

iphase1= engine.indexPhase("Kaolinite_long_phase_name");
print(iphase1, " ", engine.phaseName(iphase1))
iphase2= engine.indexPhase("Kaolinite");
print(iphase2, " ", engine.phaseName(iphase2))

#indexSpecies dependen on phase optional
print(engine.speciesAmount(ispecies1), " ", engine.speciesAmount(ispecies2))
engine.setSpeciesAmount("Anorthite_long_substance_name", 10)
engine.setSpeciesAmount("Anorthite_long_substance_name", 20, "Plagioclase")
print(engine.speciesAmount(ispecies1), " ", engine.speciesAmount(ispecies2))

print("Albite indexes\n", engine.indexSpeciesMap("Albite"))
print("Quartz indexes\n", engine.indexSpeciesMap("Quartz"))

# warning messages required
print(engine.phaseSpecificVolume(2))
print(engine.phaseSpecificVolume(30))
print(engine.phaseSpecificVolume(-1))
