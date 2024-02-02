from xgems import *
from numpy import *
from pyGEMS import GEMS

#connect all loggers
#update_loggers(True, "test_demo1.log", 0)

#engine = GEMS2("resources/CalciteBC/CalciteBC-dat.lst")
engine = GEMS2("resources/Thermo-time-in/series1-dat.lst")
element_names = engine.element_names
species_names = engine.species_names
aq_species_names = engine.species_in_phase['aq_gen']
phase_names = engine.phase_names

#print engine properties
print("T", engine.T)
print("P", engine.P)
print("b", engine.b)

print("nelements", engine.nelements)
print("nphases", engine.nphases)
print("nspecies", engine.nspecies)

print("element_names", engine.element_names)
print("species_names", engine.species_names)
print("phase_names", engine.phase_names)
print("aq_phase_symbol", engine.aq_phase_symbol)
print("gas_phase_symbol", engine.gas_phase_symbol)
print("bulk_composition", engine.bulk_composition)
print("element_molar_masses", engine.element_molar_masses)

print("species_in_phase", engine.species_in_phase)
out = engine.species_in_phase
for name in phase_names:  #print fixed order to compare
    print(name, out[name])
print("species_charges", engine.species_charges)
out = engine.species_charges
for name in species_names:
    print(name, out[name])
print("species_molar_mass", engine.species_molar_mass)
out = engine.species_molar_mass
for name in species_names:
    print(name, out[name])
print("species_molar_volumes", engine.species_molar_volumes)
out = engine.species_molar_volumes
for name in species_names:
    print(name, out[name])

print("pH", engine.pH)
print("pE", engine.pE)
print("ionic_strength", engine.ionic_strength)
print("system_volume", engine.system_volume)
print("system_mass", engine.system_mass)
print("phases_molar_volume", engine.phases_molar_volume)
out = engine.phases_molar_volume
for name in phase_names:
    print(name, out[name])
print("phase_sat_indices", engine.phase_sat_indices)
out = engine.phase_sat_indices
for name in phase_names:
    print(name, out[name])

print("aq_elements_molarity", engine.aq_elements_molarity)
print("aq_elements_molality", engine.aq_elements_molality)
print("aq_species_molarity", engine.aq_species_molarity)
out = engine.aq_species_molarity
for name in aq_species_names:
    print(name, out[name])
print("aq_species_molality", engine.aq_species_molality)
out = engine.aq_species_molality
for name in aq_species_names[:-1]:
    print(name, out[name])
print("aq_elements_moles", engine.aq_elements_moles)
print("solids_elements_moles", engine.solids_elements_moles)
print("phases_elements_moles", engine.phases_elements_moles)
out = engine.phases_elements_moles
for name in phase_names:
    print(name, out[name])

print("phases_moles", engine.phases_moles)
out = engine.phases_moles
for name in phase_names:
    print(name, out[name])

print("species_moles", engine.species_moles)
out = engine.species_moles
for name in species_names:
    print(name, out[name])
print("species_ln_activities", engine.species_ln_activities)
out = engine.species_ln_activities
for name in species_names:
    print(name, out[name])
print("species_ln_activity_coefficients", engine.species_ln_activity_coefficients)
out = engine.species_ln_activity_coefficients
for name in species_names:
    print(name, out[name])

print("phase_species_moles", engine.phase_species_moles('aq_gen'))
print("solids_mass_frac", engine.solids_mass_frac)
out = engine.solids_mass_frac
for name in phase_names:
    print(name, out[name])
print("solids_volume_frac", engine.solids_volume_frac)
print("aq_volume_frac", engine.aq_volume_frac)
print("phases_volume", engine.phases_volume)
out = engine.phases_volume
for name in phase_names:
    print(name, out[name])
print("phases_mass", engine.phases_mass)
out = engine.phases_mass
for name in phase_names:
    print(name, out[name])
print("phases_volume_frac", engine.phases_volume_frac)
out = engine.phases_volume_frac
for name in phase_names:
    print(name, out[name])


