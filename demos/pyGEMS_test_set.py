from xgems import *
from numpy import *
from pyGEMS import GEMS

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

engine = GEMS2("resources/Thermo-time-in/series1-dat.lst")
element_names = engine.element_names
species_names = engine.species_names
phase_names = engine.phase_names
print("element_names", engine.element_names)
print("species_names", engine.species_names)
print("phase_names", engine.phase_names)

# set functions -------------------------------------------------

b_from_formula1 = engine.get_b_from_formula( {"H":2,"O":1}, 1, "kg")
print("b_from_formula", b_from_formula1)
b_from_formula2 = engine.get_b_from_formula( {"H":2,"O":1})
print("b_from_formula", b_from_formula2)

engine.clear()
print("\nbulk_composition", engine.vector_b)
engine.add_multiple_species_amt( { 'HCl@':0.01, 'H2@':2 }, "moles")
engine.add_species_amt( 'H2O@', 0.01, "kg")
engine.add_element_amt( 'Al', 0.3, "moles")
engine.add_multiple_elements_amt( { 'Na':1.013077, 'Si':1.013077 }, "moles")
engine.add_amt_from_formula( { 'K':2, 'O':1 }, 4.108*1e-3, "kg")
print("\nbulk_composition", engine.bulk_composition)

values_map = {'Al': 0.2944257, 'Ca': 1.78325e-08, 'Cl': 0.01, 'H': 1.110267, 'K': 0.08722239, 'Na': 0.1783473, 'O': 3.15066, 'Si': 1.013077, 'Zz': 0.0}

engine.set_bulk_composition(values_map)
print("equilibrate = ", engine.equilibrate() )
print("phases_moles", engine.phases_moles)
out = engine.phases_moles
for name in phase_names:
    print(name, out[name])
print("species_moles", engine.species_moles)
out = engine.species_moles
for name in species_names:
    print(name, out[name])

