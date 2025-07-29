#include <map>
#include "ChemicalEngineMaps.hpp"


namespace xGEMS {

static std::map<int, std::string> _status_encoder = {
    { 0, "No GEM re-calculation needed" },
    { 1, "Need GEM calculation with LPP (automatic) initial approximation (AIA)"},
    { 2, "OK after GEM calculation with LPP AIA"},
    { 3, "Bad (not fully trustful) result after GEM calculation with LPP AIA"},
    { 4, "Failure (no result) in GEM calculation with LPP AIA"},
    { 5, "Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation (full DATABR lists only)"},
    { 6, "OK after GEM calculation with SIA"},
    { 7, "Bad (not fully trustful) result after GEM calculation with SIA"},
    { 8, "Failure (no result) in GEM calculation with SIA"},
    { 9, "Terminal error has occurred in GEMS3K (e.g. memory corruption). Restart is required."},
};


ChemicalEngineMaps::ChemicalEngineMaps(const std::string &inputfile, bool reset_calc, bool coldstart):
    input_file(inputfile), gem(inputfile)
{
    T = gem.temperature();
    P = gem.pressure();
    b_amounts = gem.elementAmounts();

    if( coldstart ) {
        cold_start();
    }

    equilibrate();

    m_aq_phase_symbol = gem.aqueousPhaseName();
    m_gas_phase_symbol = gem.gasPhaseName();

    auto elemolarmass = gem.elementMolarMasses();
    for(Index i = 0; i < nelements(); ++i) {
        m_element_names.push_back(gem.elementName(i));
        m_element_molar_masses[gem.elementName(i)] = elemolarmass[i];
    }

    auto molar_mass = gem.speciesMolarMasses();
    for(Index i = 0; i < nspecies(); ++i) {
        auto specname = gem.speciesName(i);
        m_species_names.push_back(specname);
        m_species_charges[specname] = gem.speciesCharge(i);
        m_species_molar_mass[specname] = molar_mass[i];
        m_species_molar_volumes[specname] = gem.standardMolarVolume(i);
    }

    for(Index i = 0; i < nphases(); ++i) {
        m_phase_names.push_back(gem.phaseName(i));
    }
    for(Index i = 0; i < nphases(); ++i) {
        if( gem.numSpeciesInPhase(i) > 0) {
            m_species_in_phase[m_phase_names[i]] = std::vector( m_species_names.begin()+gem.indexFirstSpeciesInPhase(i),
                                                                m_species_names.begin()+gem.indexFirstSpeciesInPhase(i)+gem.numSpeciesInPhase(i));
        }
    }

    //auto formulaMatrix = gem.formulaMatrix().T();
    if( reset_calc ) {
        clear();
    }
}

auto ChemicalEngineMaps::equilibrate() -> std::string
{
    auto outcode= gem.equilibrate(T,P,b_amounts);
    return _status_encoder[outcode];
}

auto ChemicalEngineMaps::clear_vector(Vector& bb, double min_amount) -> void
{
    if( min_amount > 0 ) {
        for(Index i = 0; i < nelements(); ++i) {
            if( m_element_names[i] == "Zz") {
                bb[i] = 0.0;
            }
            else {
                bb[i] = min_amount;
            }
        }
    }
}

auto ChemicalEngineMaps::clear(double min_amount) -> void
{
    clear_vector(b_amounts, min_amount);
}


auto ChemicalEngineMaps::set_species_G0(std::string symbol, double value) -> void
{
    gem.setStandardMolarGibbsEnergy(symbol, value);
}

// return input bulk elemental composition (vector b) in moles
auto ChemicalEngineMaps::bulk_composition() -> ValuesMap
{
    return to_map( m_element_names, b_amounts );
}

// aq solution composition in mol/L aq solution
auto ChemicalEngineMaps::aq_elements_molarity() -> ValuesMap
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(m_aq_phase_symbol);
    if(aq_index < nphases() ) {
        auto moles_elements = gem.elementAmountsInPhase(aq_index);
        for(Index i = 0; i < nelements(); ++i) {
            out[gem.elementName(i)] = moles_elements[i] / (gem.phaseVolume(aq_index)*1000);// volume from m3 to L
        }
    }
    return out;
}


// aq solution elemental composition in mol/kgH2O
auto ChemicalEngineMaps::aq_elements_molality() -> ValuesMap
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(m_aq_phase_symbol);
    if(aq_index < nphases() ) {
        auto H2Oindex = gem.numSpeciesInPhase(aq_index)-1;
        auto H2Oamount = gem.speciesAmounts()[H2Oindex];
        auto H2Ommass = gem.speciesMolarMasses()[H2Oindex];
        auto H2Omass = H2Oamount*H2Ommass/1000; // in kg
        auto moles_elements = gem.elementAmountsInPhase(aq_index);
        for(Index i = 0; i < nelements(); ++i) {
            out[gem.elementName(i)] = moles_elements[i] / H2Omass;
        }
    }
    return out;
}

// aq solution composition in mol/L of aqueous species
auto ChemicalEngineMaps::aq_species_molarity() -> ValuesMap
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(m_aq_phase_symbol);
    if(aq_index < nphases() ) {
        auto moles_species = gem.speciesAmounts();
        for(Index i = 0; i < nspecies(); ++i) {
            out[gem.speciesName(i)] =  moles_species[i] / (gem.phaseVolume(aq_index)*1000); // volume from m3 to L
        }
    }
    return out;
}

// aq solution composition in mol/kg H2O of aqueous species (speciation)
auto ChemicalEngineMaps::aq_species_molality() -> ValuesMap
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(m_aq_phase_symbol);
    if(aq_index < nphases() ) {
        auto molalities =  gem.speciesMolalities();
        for(Index i = gem.indexFirstSpeciesInPhase(aq_index);
            i < gem.indexFirstSpeciesInPhase(aq_index)+gem.numSpeciesInPhase(aq_index)-1; ++i) {
            out[gem.speciesName(i)] =  molalities[i];
        }
    }
    return out;
}


// aq solution elements amount in moles
auto ChemicalEngineMaps::aq_elements_moles() -> ValuesMap
{
    auto aq_index = gem.indexPhase(m_aq_phase_symbol);
    if(aq_index < nphases() ) {
        return to_map( m_element_names,  gem.elementAmountsInPhase(aq_index) );
    }
    return {};
}

// set input bulk elemental composition (vector b) in moles
auto ChemicalEngineMaps::set_bulk_composition(ValuesMap b_input, double min_amount) -> void
{
    for(Index i = 0; i < nelements(); ++i) {
        if( b_input.find(m_element_names[i]) != b_input.end() ) {
            b_amounts[i] = b_input[m_element_names[i]];
        }
        else if(b_amounts[i] < min_amount) {
            b_amounts[i] = min_amount;
        }
        else if(m_element_names[i] == "Zz") {
            b_amounts[i] = 0.0;
        }

    }
}

//  Removes bulk elemental aqueous solution composition from vector b
//  be careful as this will also remove water i.e H+ and OH-
//  Not quite clear what this access method really does (DK)
auto ChemicalEngineMaps::reset_aq_composition(double min_amount) -> void
{
    auto peamt = gem.phaseAmounts();
    auto aq_index = gem.indexPhase(m_aq_phase_symbol);
    if(aq_index < nphases() ) {
        if( peamt[aq_index] > min_amount ) {
            auto b_aqup = gem.elementAmountsInPhase(aq_index);
            b_amounts -= b_aqup;
        }
    }
    for(Index i = 0; i < nelements(); ++i) {
        if(b_amounts[i] < min_amount) {
            b_amounts[i] = min_amount;
        }
    }
}


// return a dictionary containing mole amounts of elements in all solids together
auto ChemicalEngineMaps::solids_elements_moles(double min_amount_phase, double min_amount_element) -> ValuesMap
{
    Vector  b_solid = gem.elementAmounts();
    auto peamt = gem.phaseAmounts();
    auto aqupx = gem.indexPhase(m_aq_phase_symbol);
    if(aqupx < nphases() ) {

        if( peamt[aqupx] > min_amount_phase) {
            auto b_aqup = gem.elementAmountsInPhase(aqupx);
            b_solid -= b_aqup;
        }
    }
    auto gaspx = gem.indexPhase(m_gas_phase_symbol);
    if( gaspx < nphases() ) {
        if(peamt[gaspx] > min_amount_phase) {
            auto b_gasp = gem.elementAmountsInPhase(gaspx);
            b_solid -= b_gasp;
        }
    }
    ValuesMap out;
    for(Index i = 0; i < nelements(); ++i) {
        out[m_element_names[i]] = (b_solid[i] < min_amount_element ? 0.0 : b_solid[i]);
    }
    return out;
}


// return a dictionary (table) containing amounts of elements in phases in moles
auto ChemicalEngineMaps::phases_elements_moles() -> PhaseValuesMap
{
    PhaseValuesMap out;
    for(Index k = 0; k < nphases(); ++k) {
        auto peamt = gem.elementAmountsInPhase(k);
        ValuesMap dictelems;
        for(Index i = 0; i < nelements(); ++i) {
            dictelems[m_element_names[i]] = peamt[i];
        }
        out[m_phase_names[k]] = dictelems;
    }
    return out;
}


// return phases amounts in moles
auto ChemicalEngineMaps::phases_moles() -> ValuesMap
{
    return to_map( m_phase_names,  gem.phaseAmounts() );
}

// returns all species amounts in moles
auto ChemicalEngineMaps::species_moles() -> ValuesMap
{
    return to_map( m_species_names,  gem.speciesAmounts() );
}

// returns species ln(activities)
auto ChemicalEngineMaps::species_ln_activities() -> ValuesMap
{
    return to_map( m_species_names,  gem.lnActivities() );
}

// returns species ln(activity_coefficient)
auto ChemicalEngineMaps::species_ln_activity_coefficients() -> ValuesMap
{
    return to_map( m_species_names,  gem.lnActivityCoefficients() );
}

// returns the upper limits for the species
auto ChemicalEngineMaps::species_upper_bounds() -> ValuesMap
{
    return to_map( m_species_names,  gem.speciesUpperLimits() );
}

// returns the lower limits for the species
auto ChemicalEngineMaps::species_lower_bounds() -> ValuesMap
{
    return to_map( m_species_names,  gem.speciesLowerLimits() );
}


// returns all species amounts in moles
auto ChemicalEngineMaps::phase_species_moles() -> PhaseValuesMap
{
    return to_phase_species_map( gem.speciesAmounts() );
}

// returns species ln(activities)
auto ChemicalEngineMaps::phase_species_ln_activities() -> PhaseValuesMap
{
    return to_phase_species_map( gem.lnActivities() );
}

// returns species ln(activity_coefficient)
auto ChemicalEngineMaps::phase_species_ln_activity_coefficients() -> PhaseValuesMap
{
    return to_phase_species_map( gem.lnActivityCoefficients() );
}

// returns the upper limits for the species
auto ChemicalEngineMaps::phase_species_upper_bounds() -> PhaseValuesMap
{
    return to_phase_species_map( gem.speciesUpperLimits() );
}

// returns the lower limits for the species
auto ChemicalEngineMaps::phase_species_lower_bounds() -> PhaseValuesMap
{
    return to_phase_species_map( gem.speciesLowerLimits() );
}

auto ChemicalEngineMaps::to_map(const std::vector<std::string> &names, Vector values) -> ValuesMap
{
    ValuesMap out;
    for(size_t j = 0; j < names.size(); ++j) {
        out[names[j]]=values[j];
    }
    return out;
}

auto ChemicalEngineMaps::to_phase_species_map(Vector values)  -> PhaseValuesMap
{
    PhaseValuesMap phase_out;
    Index jk = 0, j = 0;
    for(Index k = 0; k < nphases(); ++k) {
        ValuesMap out;
        for( jk = 0; jk < gem.numSpeciesInPhase(k); ++j, ++jk) {
            out[m_species_names[j]]=values[j];
        }
        phase_out[m_phase_names[k]] = out;
    }
    return phase_out;
}

// returns species in phase in moles
auto ChemicalEngineMaps::phase_species_moles(std::string phase_symbol) -> ValuesMap
{
    ValuesMap out;
    auto index = gem.indexPhase(phase_symbol);
    if( index < nphases() ) {
        auto amounts =  gem.speciesAmounts();
        for(Index i = gem.indexFirstSpeciesInPhase(index);
            i < gem.indexFirstSpeciesInPhase(index)+gem.numSpeciesInPhase(index); ++i) {
            out[m_species_names[i]] =  amounts[i];
        }
    }
    return out;
}


// mass(phase)/mass(system) ratios for [solid] phases
auto ChemicalEngineMaps::solids_mass_frac() -> ValuesMap
{
    Vector mfrac = gem.phaseMasses();
    auto sum = mfrac.sum();
    mfrac = mfrac/sum;
    return to_map( m_phase_names, mfrac );
}

// volume(phase)/volume(total) ratio for solid phases
auto ChemicalEngineMaps::solids_volume_frac() -> ValuesMap
{
    auto out = phases_volume_frac();
    out.erase(m_aq_phase_symbol);
    out.erase(m_gas_phase_symbol);
    return out;
}

// Volume fraction of aqueous phase in the system
auto ChemicalEngineMaps::aq_volume_frac() -> double
{
    auto out = phases_volume_frac();
    return out[m_aq_phase_symbol];
}

// returns a dict. with phases and their absolute volume in m3
auto ChemicalEngineMaps::phases_volume() -> ValuesMap
{
    return to_map( m_phase_names, gem.phaseVolumes() );
}

// returns a dict. with phases and their mass in kg
auto ChemicalEngineMaps::phases_mass() -> ValuesMap
{
    return to_map( m_phase_names, gem.phaseMasses() );
}

// returns a dict. with phases and their volume fractions in the system
auto ChemicalEngineMaps::phases_volume_frac() -> ValuesMap
{
    auto vfrac = gem.phaseVolumes()/system_volume();
    return to_map( m_phase_names, vfrac );
}

//  add species amount in the system useful for adding aqueous solution composition
//  units= moles, kg, m3
auto ChemicalEngineMaps::add_multiple_species_amt(const ValuesMap& input_dict, const std::string& units) -> void
{
    for(const auto& el: input_dict) {
        add_species_amt(el.first, el.second, units);
    }
}

// add species amount in the system useful for adding aqueous solution composition
// units= moles, kg, m3
auto ChemicalEngineMaps::add_species_amt(const std::string& species, double val, const std::string& units) -> void
{
    auto species_idx =gem.indexSpecies(species);
    if( units == "kg") {
        val/=m_species_molar_mass[species];
    }
    if( units == "m3") {
        val/=m_species_molar_volumes[species];
    }
    MatrixConstRef W = gem.formulaMatrix();
    MatrixConstRef Wp = W(all, species_idx);
    b_amounts += Wp*val;
}

// add element amount in the system units = moles, kg
auto ChemicalEngineMaps::add_element_amt(const std::string& element_name, double val, const std::string& units) -> void
{
    if( units  == "kg" ) {
        val /= m_element_molar_masses[element_name];
    }
    auto el_index = gem.indexElement(element_name);
    b_amounts[el_index] += val;
}

//  add elements amount in the system useful for adding aqueous solution composition
//  units= moles,kg
auto ChemicalEngineMaps::add_multiple_elements_amt(const ValuesMap& input_dict, const std::string& units) -> void
{
    for(const auto& el: input_dict) {
        add_element_amt(el.first, el.second, units);
    }
}

// add element amount using user defined formula, units = moles,kg
auto ChemicalEngineMaps::add_amt_from_formula(const ValuesMap& formula, double val, const std::string& units) -> void
{
    if( units  == "kg" ) {
        double molarmass =0.0;
        for( const auto& element: formula ) {
            molarmass += element.second * m_element_molar_masses[element.first];
        }
        val/=molarmass;
    }
    for( const auto& element: formula ) {
        add_element_amt(element.first, val * element.second);
    }
}


// returns a bulk vector b from user-defined formula (as dict. {"H":2,"O":1} )
// and amount of the formula [object] in units of 'moles' or 'kg'
auto ChemicalEngineMaps::get_b_from_formula(const ValuesMap& formula, double val, const std::string& units, double min_amount) -> Vector
{
    Vector bx(b_amounts.size()); //   bx = [v for v in self.b]
    clear_vector(bx, min_amount);

    if( units  == "kg" ) {
        double molarmass =0.0;
        for( const auto& element: formula ) {
            molarmass += element.second * m_element_molar_masses[element.first];
        }
        val/=molarmass;
    }
    for( const auto& element: formula ) {
        auto el_index = gem.indexElement(element.first);
        bx[el_index] += val * element.second;
    }
    return bx;
}

//     constrain species amount to a specified lower bound, units= moles,kg,m3
auto ChemicalEngineMaps::set_multiple_species_lower_bound(const ValuesMap& input_dict, const std::string& units) -> void
{
    for(const auto& el: input_dict) {
        set_species_lower_bound(el.first, el.second, units);
    }
}

//     constrain species amount to a specified lower bound, units= moles,kg,m3
auto ChemicalEngineMaps::set_multiple_species_upper_bound(const ValuesMap& input_dict, const std::string& units) -> void
{
    for(const auto& el: input_dict) {
        set_species_upper_bound(el.first, el.second, units);
    }
}

//  constrain species amount to a specified lower bound, units= moles,kg,m3
auto ChemicalEngineMaps::set_species_lower_bound(const std::string& species, double val, const std::string& units) -> void
{
    if( units == "kg") {
        val/=m_species_molar_mass[species];
    }
    if( units == "m3") {
        val/=m_species_molar_volumes[species];
    }
    gem.setSpeciesLowerLimit(species,val);
}

//  constrain species amount to a specified upper bound, units= moles,kg,m3
auto ChemicalEngineMaps::set_species_upper_bound(const std::string& species, double val, const std::string& units) -> void
{
    if( units == "kg") {
        val/=m_species_molar_mass[species];
    }
    if( units == "m3") {
        val/=m_species_molar_volumes[species];
    }
    gem.setSpeciesUpperLimit(species,val);
}

//  constrain species amount to a specified lower bound, units= moles,kg,m3
auto ChemicalEngineMaps::set_species_lower_bound(Index ispecies, double val, const std::string& units) -> void
{
    if( units == "kg") {
        val/=m_species_molar_mass[m_species_names[ispecies]];
    }
    if( units == "m3") {
        val/=m_species_molar_volumes[m_species_names[ispecies]];
    }
    gem.setSpeciesLowerLimit(ispecies,val);
}

//  constrain species amount to a specified upper bound, units= moles,kg,m3
auto ChemicalEngineMaps::set_species_upper_bound(Index ispecies, double val, const std::string& units) -> void
{
    if( units == "kg") {
        val/=m_species_molar_mass[m_species_names[ispecies]];
    }
    if( units == "m3") {
        val/=m_species_molar_volumes[m_species_names[ispecies]];
    }
    gem.setSpeciesUpperLimit(ispecies,val);
}


// supresses a phase in GEM calculation
auto ChemicalEngineMaps::supress_phase(const std::string& phase_name, double min_amount, double max_amount) -> void
{
    supress_multiple_species(m_species_in_phase[phase_name], min_amount, max_amount);
}

// supresses multiple phase in calculation as given in phase names list
auto ChemicalEngineMaps::supress_multiple_phases(const std::vector<std::string>& phase_name_list, double min_amount, double max_amount) -> void
{
    for( const auto& phase: phase_name_list) {
        supress_phase(phase, min_amount, max_amount);
    }
}

// supresses species in calculation
auto ChemicalEngineMaps::supress_species(const std::string& species_name, double min_amount, double max_amount) -> void
{
    set_species_lower_bound(species_name, min_amount);
    set_species_upper_bound(species_name, max_amount);
}

// supresses multiple species in GEM calculation as given in species name list
auto ChemicalEngineMaps::supress_multiple_species(const std::vector<std::string>& species_list, double min_amount, double max_amount) -> void
{
    for( const auto& species: species_list) {
        supress_species(species, min_amount, max_amount);
    }
}

// activate supressed phase
auto ChemicalEngineMaps::activate_phase(const std::string& phase_name) -> void
{
    activate_multiple_species(m_species_in_phase[phase_name]);
}

// activate multiple supressed phases given in list
auto ChemicalEngineMaps::activate_multiple_phases(const std::vector<std::string>& phase_name_list) -> void
{
    for( const auto& phase: phase_name_list) {
        activate_phase(phase);
    }
}

// activate multiple supressed species given in the list
auto ChemicalEngineMaps::activate_multiple_species(const std::vector<std::string>& species_list) -> void
{
    for( const auto& species: species_list) {
        activate_species(species);
    }
}

// activate a supressed species in phase
auto ChemicalEngineMaps::activate_species(const std::string& species_name) -> void
{
    set_species_lower_bound(species_name, 0);
    set_species_upper_bound(species_name, 1e6);
}


// returns pH of the solution
auto ChemicalEngineMaps::pH() -> double
{
    return gem.pH();
}

// returns pE of the solution
auto ChemicalEngineMaps::pE() -> double
{
    return gem.pe();
}

// returns ionic strength of the solution
auto ChemicalEngineMaps::ionic_strength() -> double
{
    return gem.ionicStrength();
}

// returns volume of the system in m3
auto ChemicalEngineMaps::system_volume() -> double
{
    return gem.systemVolume();
}

// returns mass of the system in kg
auto ChemicalEngineMaps::system_mass() -> double
{
    return gem.systemMass();
}

// returns molar volume of phases in m3/mol
auto ChemicalEngineMaps::phases_molar_volume() -> ValuesMap
{
    ValuesMap phase_mvol;
    for(Index i = 0; i < nphases(); ++i) {
        phase_mvol[m_phase_names[i]] = gem.phaseMolarVolume(i);
    }
    return phase_mvol;
}

// returns saturation indices of phases
auto ChemicalEngineMaps::phase_sat_indices() -> ValuesMap
{
    return to_map(m_phase_names, gem.phaseSatIndices());
}

}
