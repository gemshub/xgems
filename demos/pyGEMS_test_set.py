from xgems import *
from numpy import *

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

engine = GEMS2("resources/Thermo-time-in/series1-dat.lst")
element_names = engine.element_names
species_names = engine.species_names
phase_names = engine.phase_names
print("element_names\n", engine.element_names)
print("species_names\n", engine.species_names)
print("phase_names\n", engine.phase_names)

# set functions -------------------------------------------------

b_from_formula1 = engine.get_b_from_formula( {"H":2,"O":1}, 1, "kg")
print("\nb_from_formula", b_from_formula1)
b_from_formula2 = engine.get_b_from_formula( {"H":2,"O":1})
print("\nb_from_formula", b_from_formula2)

engine.clear()
print("\nempty bulk_composition\n", engine.vector_b)
engine.add_multiple_species_amt( { 'HCl@':0.01, 'H2@':2 }, "moles")
engine.add_species_amt( 'H2O@', 0.01, "kg")
engine.add_element_amt( 'Al', 0.3, "moles")
engine.add_multiple_elements_amt( { 'Na':1.013077, 'Si':1.013077 }, "moles")
engine.add_amt_from_formula( { 'K':2, 'O':1 }, 4.108*1e-3, "kg")
print("\nbulk_composition\n", engine.bulk_composition)


values_map = {'Al': 0.2944257, 'Ca': 1.78325e-08, 'Cl': 0.01, 'H': 1.110267, 'K': 0.08722239, 'Na': 0.1783473, 'O': 3.15066, 'Si': 1.013077, 'Zz': 0.0}
engine.set_bulk_composition(values_map)
print("equilibrate = ", engine.equilibrate() )

out1 = engine.phase_sat_indices
out2 = engine.phases_molar_volume
out3 = engine.phases_moles
out4 = engine.phases_mass
out5 = engine.phases_volume
out6 = engine.phases_volume_frac
out7 = engine.solids_mass_frac

print("\n           phase; phase_sat_indices; phases_molar_volume;  phases_moles;     phases_mass;   phases_volume; phases_volume_frac; solids_mass_frac;")
for name in phase_names:
    print("{:16s}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g}; {:16.8g};".format(
           name, out1[name], out2[name], out3[name], out4[name], out5[name], out6[name], out7[name]))

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


