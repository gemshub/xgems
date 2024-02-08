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
print("aq_phase_symbol", engine.aq_phase_symbol)
print("gas_phase_symbol", engine.gas_phase_symbol)
print("pH", engine.pH)
print("pE", engine.pE)
print("ionic_strength", engine.ionic_strength)
print("system_volume", engine.system_volume)
print("system_mass", engine.system_mass)
print("aq_volume_frac", engine.aq_volume_frac)

#elements
print("\nelements---------------\n")
print("nelements", engine.nelements)
print("element_names", engine.element_names)
print("bulk_composition", engine.bulk_composition)
print("element_molar_masses", engine.element_molar_masses)

#phases
print("\nphases---------------\n")
print("nphases", engine.nphases)
print("phase_names", engine.phase_names)
print("species_in_phase", engine.species_in_phase)
out = engine.species_in_phase
for name in phase_names:  #print fixed order to compare
    print(name, out[name])

print("\nphase_sat_indices", engine.phase_sat_indices)
print("phases_molar_volume", engine.phases_molar_volume)
print("phases_moles", engine.phases_moles)
print("phases_mass", engine.phases_mass)
print("phases_volume", engine.phases_volume)
print("phases_volume_frac", engine.phases_volume_frac)

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

print("\nphases_elements_moles", engine.phases_elements_moles)
out = engine.phases_elements_moles
for name in phase_names:
    print(name, out[name])

print("\nphase_species_moles ('aq_gen')", engine.phase_species_moles('aq_gen'))

print("\nsolids_elements_moles", engine.solids_elements_moles)
print("solids_volume_frac", engine.solids_volume_frac)
print("solids_mass_frac", engine.solids_mass_frac)

#species
print("\nspecies---------------\n")
print("nspecies", engine.nspecies)
print("species_names", engine.species_names)

print("species_charges", engine.species_charges)
print("species_molar_mass", engine.species_molar_mass)
print("species_molar_volumes", engine.species_molar_volumes)
print("species_moles", engine.species_moles)
print("species_ln_activities", engine.species_ln_activities)
print("species_ln_activity_coefficients", engine.species_ln_activity_coefficients)

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

#aq phase
print("\naq phase---------------\n")
print("aq_elements_moles", engine.aq_elements_moles)
print("aq_elements_molarity", engine.aq_elements_molarity)
print("aq_elements_molality", engine.aq_elements_molality)
print("aq_species_molarity", engine.aq_species_molarity)
print("aq_species_molality", engine.aq_species_molality)

out1 = engine.aq_species_molarity
out2 = engine.aq_species_molality

print("\n     specie; aq_species_molarity;  aq_species_molality;")
for name in aq_species_names[:-1]:
    print("{:12s}; {:16.8g}; {:16.8g}; ".format( name, out1[name], out2[name]))






