#ifndef GEMSENGINE_H
#define GEMSENGINE_H

#include <vector>
#include <map>
#include "ChemicalEngine.hpp"

namespace xGEMS {

using ValuesMap = std::map<std::string, double>;

/// Gems interface in calculator format for easy  using dictionaries
class GEMSEngine
{
public:

    ///  Initialization of the calculator
    ///  @param   reset_calc: true will reset b vector to zero
    GEMSEngine(const std::string& input_file, bool reset_calc=true, bool cold_start=true);

    /// runs GEM equilibriation of the current (internally set) system
    std::string equilibrate();

    void cold_start()
    {
        gem.setColdStart();
    }

    void warm_start()
    {
        gem.setWarmStart();
    }

    // number of elements in the system
    int nelements()
    {
        return gem.numElements();
    }

    // number of phases in the system
    int nphases()
    {
        return gem.numPhases();
    }

    int nspecies()
    {
        return gem.numSpecies();
    }

    void clear(double cvalue=1e-15);

    void set_species_G0(std::string symbol, double value);

    ///  return input bulk elemental composition (vector b) in moles
    ValuesMap bulk_composition();

    /// aq solution composition in mol/L aq solution
    ValuesMap aq_elements_molarity();
    /// aq solution elemental composition in mol/kgH2O
    ValuesMap aq_elements_molality();
    /// aq solution composition in mol/L of aqueous species
    ValuesMap aq_species_molarity();
    /// aq solution composition in mol/kg H2O of aqueous species (speciation)
    ValuesMap aq_species_molality();
    /// aq solution elements amount in moles
    ValuesMap aq_elements_moles();
    /// set input bulk elemental composition (vector b) in moles
    void set_bulk_composition(ValuesMap b_input);
    ///  Removes bulk elemental aqueous solution composition from vector b
    ///  be careful as this will also remove water i.e H+ and OH-
    void reset_aq_composition();
    /// return a dictionary containing mole amounts of elements in all solids together
    ValuesMap solids_elements_moles();
    /// return a dictionary (table) containing amounts of elements in phases in moles
    std::map<std::string, ValuesMap> phases_elements_moles();

    /// return phases amounts in moles
    ValuesMap phases_moles();
    /// returns all species amounts in moles
    ValuesMap species_moles();
    /// returns species ln(activities)
    ValuesMap species_ln_activities();
    /// returns species ln(activity_coefficient)
    ValuesMap species_ln_activity_coefficients();
    /// returns species in phase in moles
    ValuesMap phase_species_moles(const std::string &phase_symbol);
    /// mass(phase)/mass(system) ratios for [solid] phases
    ValuesMap solids_mass_frac();
    /// volume(phase)/volume(total) ratio for solid phases
    ValuesMap solids_volume_frac();
    /// Volume fraction of aqueous phase in the system
    double aq_volume_frac();
    /// returns a dict. with phases and their absolute volume in m3
    ValuesMap phases_volume();
    /// returns a dict. with phases and their mass in kg
    ValuesMap phases_mass();
    /// returns a dict. with phases and their volume fractions in the system
    ValuesMap phases_volume_frac();
    ///  add species amount in the system useful for adding aqueous solution composition
    ///  units= moles, kg, m3
    void add_multiple_species_amt(const ValuesMap &input_dict, const std::string& units = "moles");

    /// add species amount in the system useful for adding aqueous solution composition
    /// units= moles, kg, m3
    void add_species_amt(const std::string &species, double val, const std::string &units = "moles");
    /// add element amount in the system units = moles, kg
    void add_element_amt(const std::string &element_name, double val, const std::string &units = "moles");
    ///  add elements amount in the system useful for adding aqueous solution composition
    ///  units= moles,kg
    void add_multiple_elements_amt(const ValuesMap &input_dict, const std::string &units = "moles");
    /// add element amount using user defined formula, units = moles,kg
    void add_amt_from_formula(const ValuesMap &formula, double val, const std::string &units = "moles");
    /// returns a bulk vector b from user-defined formula (as dict. {"H":2,"O":1} )
    /// and amount of the formula [object] in units of 'moles' or 'kg'
    Vector get_b_from_formula(const ValuesMap &formula, double val = 1, const std::string &units = "moles");

    ///  constrain species amount to a specified lower bound, units= moles,kg,m3
    void set_multiple_species_lower_bound(const ValuesMap &input_dict, const std::string &units = "moles");
    ///  constrain species amount to a specified lower bound, units= moles,kg,m3
    void set_multiple_species_upper_bound(const ValuesMap &input_dict, const std::string &units = "moles");
    ///  constrain species amount to a specified lower bound, units= moles,kg,m3
    void set_species_lower_bound(const std::string &species, double val, const std::string &units= "moles");
    ///  constrain species amount to a specified upper bound, units= moles,kg,m3
    void set_species_upper_bound(const std::string &species, double val, const std::string &units= "moles");
    /// supresses a phase in GEM calculation
    void supress_phase(const std::string &phase_name);
    /// supresses multiple phase in calculation as given in phase names list
    void supress_multiple_phases(const std::vector<std::string> &phase_name_list);
    /// supresses species in calculation
    void supress_species(const std::string &species_name);
    /// supresses multiple species in GEM calculation as given in species name list
    void supress_multiple_species(const std::vector<std::string> &species_list);
    /// activate supressed phase
    void activate_phase(const std::string &phase_name);
    /// activate multiple supressed phases given in list
    void activate_multiple_phases(const std::vector<std::string> &phase_name_list);
    /// activate multiple supressed species given in the list
    void activate_multiple_species(const std::vector<std::string> &species_list);
    /// activate a supressed species in phase
    void activate_species(const std::string &species_name);

    /// returns pH of the solution
    double pH();
    /// returns pE of the solution
    double pE();
    /// returns ionic strength of the solution
    double ionic_strength();
    /// returns volume of the system in m3
    double system_volume();
    /// returns mass of the system in kg
    double system_mass();
    /// returns molar volume of phases in m3/mol
    ValuesMap phases_molar_volume();
    /// returns saturation indices of phases
    ValuesMap phase_sat_indices();

protected:

    std::string input_file;
    ChemicalEngine gem;
    double T; // K
    double P; // Pa
    Vector b_amounts; // moles

    std::vector<std::string> element_names;
    std::vector<std::string> species_names;
    std::vector<std::string> phase_names;

    std::string aq_phase_symbol;
    std::string gas_phase_symbol;

    /// dictionary containing elements and their molar masses
    ValuesMap element_molar_masses;
    /// dictionary containing species in phase
    std::map<std::string, std::vector<std::string>> species_in_phase;
    ///
    ValuesMap species_charges;
    /// dictionary containing species and their molar masses
    ValuesMap species_molar_mass;
    ///  dictionary containing species and their molar volumes
    ValuesMap species_molar_volumes;



    ValuesMap to_map( const std::vector<std::string>& names,  Vector values )
    {
        ValuesMap out;
        for(Index j = 0; j < names.size(); ++j) {
            out[names[j]]=values[j];
        }
        return out;
    }


};

}
#endif // GEMSENGINE_H
