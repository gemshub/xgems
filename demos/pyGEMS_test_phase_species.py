from xgems import *
from numpy import *

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

engine = GEMS2("resources/Thermo-time-in/series1-dat.lst")
phase_names = engine.phase_names

print("phase_names\n", engine.phase_names)
print("\nspecies_in_phase\n", engine.species_in_phase)
out = engine.species_in_phase
for name in phase_names:
    print(name, out[name])

print("\nphase_species_upper_bounds\n", engine.phase_species_upper_bounds)
print("\nphase_species_lower_bounds\n", engine.phase_species_lower_bounds)
print("\nphase_species_moles")
out = engine.phase_species_moles()
for name in phase_names:
    print(name, out[name])

print("\nphase_species_ln_activities\n", engine.phase_species_ln_activities)

print("\nphase_species_ln_activity_coefficients")
out = engine.phase_species_ln_activity_coefficients
for name in phase_names:
    print(name, out[name])

