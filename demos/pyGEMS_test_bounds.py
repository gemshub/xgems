from xgems import *
from numpy import *

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

engine = GEMS2("resources/CalciteBC/CalciteBC-dat.lst")
element_names = engine.element_names
species_names = engine.species_names
phase_names = engine.phase_names

print("test set bounds---------------------------------------")

print("species_bounds")
out1 = engine.species_lower_bounds
out2 = engine.species_upper_bounds
for name in species_names:
    print(name, out1[name], out2[name])


engine.set_species_lower_bound( 8, 400, "moles")
engine.set_species_upper_bound( 8, 900, "kg")
engine.set_species_lower_bound( 'Ca(HCO3)+', 200, "moles")
engine.set_species_upper_bound( 'CaOH+', 500, "kg")
engine.set_multiple_species_lower_bound( {'Mg(CO3)@':30, 'Mg(HCO3)+':40, 'Mg+2':50})
engine.set_multiple_species_upper_bound( {'Mg(CO3)@':300, 'Mg(HCO3)+':400, 'Mg+2':500})
engine.supress_phase('gas_gen')
engine.supress_multiple_phases(['Dolomite-dis', 'Tin'])
engine.supress_species('Ca(CO3)@')
engine.supress_multiple_species(['ClO4-', 'Cl-'])

print("\nspecies_bounds_supress")
out1 = engine.species_lower_bounds
out2 = engine.species_upper_bounds
for name in species_names:
    print(name, out1[name], out2[name])

engine.activate_phase('gas_gen')
engine.activate_multiple_phases(['Dolomite-dis', 'Tin'])
engine.activate_species('Ca(CO3)@')
engine.activate_multiple_species(['Ca(HCO3)+', 'CaOH+', 'Mg(CO3)@', 'Mg(HCO3)+', 'Mg+2', 'ClO4-', 'Cl-'])

print("\nspecies_bounds_activate")
out1 = engine.species_lower_bounds
out2 = engine.species_upper_bounds
for name in species_names:
    print(name, out1[name], out2[name])

