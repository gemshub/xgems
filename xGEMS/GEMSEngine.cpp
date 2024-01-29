#include <map>
#include "GEMSEngine.h"


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


GEMSEngine::GEMSEngine(const std::string &inputfile, bool reset_calc, bool coldstart):
    input_file(inputfile), gem(inputfile)
{
    T = gem.temperature();
    P = gem.pressure();
    b_amounts = gem.elementAmounts();

    if( coldstart ) {
        cold_start();
    }

    equilibrate();

    aq_phase_symbol = gem.phaseName(0);
    //self.aq_phase_symbol = 'aq_gen'               Potential problem if no aq phase in system
    //self.gas_phase_symbol = self.gem.phaseName(1) # use index 0 if no aq phase in the system, or
    gas_phase_symbol = "gas_gen";             //!!!!!!!!!!! Potential problem if no gas phase in system
    // to be done or it has a different name

    auto elemolarmass = gem.elementMolarMasses();
    for(Index i = 0; i < nelements(); ++i) {
        element_names.push_back(gem.elementName(i));
        element_molar_masses[gem.elementName(i)] = elemolarmass[i];
    }

    auto molar_mass = gem.speciesMolarMasses();
    for(Index i = 0; i < nspecies(); ++i) {
        auto specname = gem.speciesName(i);
        species_names.push_back(specname);
        species_charges[specname] = gem.speciesCharge(i);
        species_molar_mass[specname] = molar_mass[i];
        species_molar_volumes[specname] = gem.standardMolarVolume(i);
    }

    for(Index i = 0; i < nphases(); ++i) {
        phase_names.push_back(gem.phaseName(i));
    }
    for(Index i = 0; i < nphases(); ++i) {
        if( gem.numSpeciesInPhase(i) > 0) {
            species_in_phase[phase_names[i]] = std::vector( species_names.begin()+gem.indexFirstSpeciesInPhase(i),
                                                            species_names.begin()+gem.indexFirstSpeciesInPhase(i)+gem.numSpeciesInPhase(i));
        }
    }

    //auto formulaMatrix = gem.formulaMatrix().T();
    if( reset_calc ) {
        clear();
    }
}

std::string GEMSEngine::equilibrate()
{
    auto outcode= gem.equilibrate(T,P,b_amounts);
    return _status_encoder[outcode];
}

void GEMSEngine::clear(double cvalue)
{
    if( cvalue > 0 ) {
        for(Index i = 0; i < nelements(); ++i) {
            if( element_names[i] == "Zz") {
                b_amounts[i] = 0.0;
            }
            else {
                b_amounts[i] = cvalue;
            }
        }
    }
}


void GEMSEngine::set_species_G0(std::string symbol, double value)
{
    gem.setStandardMolarGibbsEnergy(symbol, value);
}

// return input bulk elemental composition (vector b) in moles
ValuesMap GEMSEngine::bulk_composition()
{
    return to_map( element_names, b_amounts );
}

// aq solution composition in mol/L aq solution
ValuesMap GEMSEngine::aq_elements_molarity()
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(aq_phase_symbol);
    if(aq_index < nphases() ) {
        auto moles_elements = gem.elementAmountsInPhase(aq_index);
        for(Index i = 0; i < nelements(); ++i) {
            out[gem.elementName(i)] = moles_elements[i] / (gem.phaseVolume(aq_index)*1000);// volume from m3 to L
        }
    }
    return out;
}


// aq solution elemental composition in mol/kgH2O
ValuesMap GEMSEngine::aq_elements_molality()
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(aq_phase_symbol);
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
ValuesMap GEMSEngine::aq_species_molarity()
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(aq_phase_symbol);
    if(aq_index < nphases() ) {
        auto moles_species = gem.speciesAmounts();
        for(Index i = 0; i < nspecies(); ++i) {
            out[gem.speciesName(i)] =  moles_species[i] / (gem.phaseVolume(aq_index)*1000); // volume from m3 to L
        }
    }
    return out;
}

// aq solution composition in mol/kg H2O of aqueous species (speciation)
ValuesMap GEMSEngine::aq_species_molality()
{
    ValuesMap out;
    auto aq_index = gem.indexPhase(aq_phase_symbol);
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
ValuesMap GEMSEngine::aq_elements_moles()
{
    auto aq_index = gem.indexPhase(aq_phase_symbol);
    if(aq_index < nphases() ) {
        return to_map( element_names,  gem.elementAmountsInPhase(aq_index) );
    }
    return {};
}

// set input bulk elemental composition (vector b) in moles
void GEMSEngine::set_bulk_composition(ValuesMap b_input)
{
    for(Index i = 0; i < nelements(); ++i) {
        if( b_input.find(element_names[i]) != b_input.end() ) {
            b_amounts[i] = b_input[element_names[i]];
        }
        else if(b_amounts[i] < 1e-15) {
            b_amounts[i] = 1e-15;
        }
        else if(element_names[i] == "Zz") {
            b_amounts[i] = 0.0;
        }

    }
}

//  Removes bulk elemental aqueous solution composition from vector b
//  be careful as this will also remove water i.e H+ and OH-
//  Not quite clear what this access method really does (DK)
void GEMSEngine::reset_aq_composition()
{
    auto peamt = gem.phaseAmounts();
    auto aq_index = gem.indexPhase(aq_phase_symbol);
    if(aq_index < nphases() ) {
        if( peamt[aq_index] > 1e-12 ) {
            auto b_aqup = gem.elementAmountsInPhase(aq_index);
            b_amounts -= b_aqup;
        }
    }
    for(Index i = 0; i < nelements(); ++i) {
        if(b_amounts[i] < 1e-12) {
            b_amounts[i] = 1e-12;
        }
    }
}


// return a dictionary containing mole amounts of elements in all solids together
ValuesMap GEMSEngine::solids_elements_moles()
{
    Vector  b_solid = gem.elementAmounts();
    auto peamt = gem.phaseAmounts();
    auto aqupx = gem.indexPhase(aq_phase_symbol);
    if(aqupx < nphases() ) {

        if( peamt[aqupx] > 1e-12) {
            auto b_aqup = gem.elementAmountsInPhase(aqupx);
            b_solid -= b_aqup;
        }
    }
    auto gaspx = gem.indexPhase(gas_phase_symbol);
    if( gaspx < nphases() ) {
        if( peamt[gaspx] > 1e-12 ) {
            auto b_gasp = gem.elementAmountsInPhase(gaspx);
            b_solid -= b_gasp;
        }
    }
    ValuesMap out;
    for(Index i = 0; i < nelements(); ++i) {
        out[element_names[i]] = ( b_solid[i] < 1e-14 ? 0.0 : b_solid[i]);
    }
    return out;
}


// return a dictionary (table) containing amounts of elements in phases in moles
std::map<std::string, ValuesMap> GEMSEngine::phases_elements_moles()
{
    std::map<std::string, ValuesMap> out;
    for(Index k = 0; k < nphases(); ++k) {
        auto peamt = gem.elementAmountsInPhase(k);
        ValuesMap dictelems;
        for(Index i = 0; i < nelements(); ++i) {
            dictelems[element_names[i]] = peamt[i];
        }
        out[phase_names[k]] = dictelems;
    }
    return out;
}


// return phases amounts in moles
ValuesMap GEMSEngine::phases_moles()
{
    return to_map( phase_names,  gem.phaseAmounts() );
}

// returns all species amounts in moles
ValuesMap GEMSEngine::species_moles()
{
    return to_map( species_names,  gem.speciesAmounts() );
}

// returns species ln(activities)
ValuesMap GEMSEngine::species_ln_activities()
{
    return to_map( species_names,  gem.lnActivities() );
}

// returns species ln(activity_coefficient)
ValuesMap GEMSEngine::species_ln_activity_coefficients()
{
    return to_map( species_names,  gem.lnActivityCoefficients() );
}


// returns species in phase in moles
ValuesMap GEMSEngine::phase_species_moles(const std::string& phase_symbol)
{
    ValuesMap out;
    auto index = gem.indexPhase(phase_symbol);
    if( index < nphases() ) {
        auto amounts =  gem.speciesAmounts();
        for(Index i = gem.indexFirstSpeciesInPhase(index);
            i < gem.indexFirstSpeciesInPhase(index)+gem.numSpeciesInPhase(index); ++i) {
            out[species_names[i]] =  amounts[i];
        }
    }
    return out;
}


// mass(phase)/mass(system) ratios for [solid] phases
ValuesMap GEMSEngine::solids_mass_frac()
{
    Vector mfrac = gem.phaseMasses();
    auto sum = mfrac.sum();
    mfrac = mfrac/sum;
    return to_map( phase_names, mfrac );
}

// volume(phase)/volume(total) ratio for solid phases
ValuesMap GEMSEngine::solids_volume_frac()
{
    auto out = phases_volume_frac();
    out.erase(aq_phase_symbol);
    out.erase(gas_phase_symbol);
    return out;
}

// Volume fraction of aqueous phase in the system
double GEMSEngine::aq_volume_frac()
{
    auto out = phases_volume_frac();
    return out[aq_phase_symbol];
}

// returns a dict. with phases and their absolute volume in m3
ValuesMap GEMSEngine::phases_volume()
{
    return to_map( phase_names, gem.phaseVolumes() );
}

// returns a dict. with phases and their mass in kg
ValuesMap GEMSEngine::phases_mass()
{
    return to_map( phase_names, gem.phaseMasses() );
}

// returns a dict. with phases and their volume fractions in the system
ValuesMap GEMSEngine::phases_volume_frac()
{
    auto vfrac = gem.phaseVolumes()/system_volume();
    return to_map( phase_names, gem.phaseMasses() );
}

//  add species amount in the system useful for adding aqueous solution composition
//  units= moles, kg, m3
void GEMSEngine::add_multiple_species_amt(const ValuesMap& input_dict, const std::string& units)
{
    for(const auto& el: input_dict) {
       add_species_amt(el.first, el.second, units);
    }
}

// add species amount in the system useful for adding aqueous solution composition
// units= moles, kg, m3
void GEMSEngine::add_species_amt(const std::string& species, double val, const std::string& units)
{
    auto species_idx =gem.indexSpecies(species);
    if( units == "kg") {
        val/=species_molar_mass[species];
    }
    if( units == "m3") {
        val/=species_molar_volumes[species];
    }
    MatrixConstRef W = gem.formulaMatrix();
    MatrixConstRef Wp = W(all, species_idx);
    b_amounts += Wp*val;
}

// add element amount in the system units = moles, kg
void GEMSEngine::add_element_amt(const std::string& element_name, double val, const std::string& units)
{
    if( units  == "kg" ) {
        val /= element_molar_masses[element_name];
    }
    auto el_index = gem.indexElement(element_name);
    b_amounts[el_index] += val;
}

//  add elements amount in the system useful for adding aqueous solution composition
//  units= moles,kg
void GEMSEngine::add_multiple_elements_amt(const ValuesMap& input_dict, const std::string& units)
{
    for(const auto& el: input_dict) {
       add_element_amt(el.first, el.second, units);
    }
}

// add element amount using user defined formula, units = moles,kg
void GEMSEngine::add_amt_from_formula(const ValuesMap& formula, double val, const std::string& units)
{
    if( units  == "kg" ) {
        double molarmass =0.0;
        for( const auto& element: formula ) {
            molarmass += element.second * element_molar_masses[element.first];
        }
        val/=molarmass;
    }
    for( const auto& element: formula ) {
        add_element_amt(element.first, val * element.second);
    }
}


// returns a bulk vector b from user-defined formula (as dict. {"H":2,"O":1} )
// and amount of the formula [object] in units of 'moles' or 'kg'
Vector GEMSEngine::get_b_from_formula(const ValuesMap& formula, double val, const std::string& units)
{
    Vector bx(b_amounts.size()); //   bx = [v for v in self.b]
    if( units  == "kg" ) {
        double molarmass =0.0;
        for( const auto& element: formula ) {
            molarmass += element.second * element_molar_masses[element.first];
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
void GEMSEngine::set_multiple_species_lower_bound(const ValuesMap& input_dict, const std::string& units)
{
    for(const auto& el: input_dict) {
       set_species_lower_bound(el.first, el.second, units);
    }
}

//     constrain species amount to a specified lower bound, units= moles,kg,m3
void GEMSEngine::set_multiple_species_upper_bound(const ValuesMap& input_dict, const std::string& units)
{
    for(const auto& el: input_dict) {
       set_species_upper_bound(el.first, el.second, units);
    }
}

//  constrain species amount to a specified lower bound, units= moles,kg,m3
void GEMSEngine::set_species_lower_bound(const std::string& species, double val, const std::string& units)
{
    auto species_idx =gem.indexSpecies(species);
    if( units == "kg") {
        val/=species_molar_mass[species];
    }
    if( units == "m3") {
        val/=species_molar_volumes[species];
    }
    gem.setSpeciesLowerLimit(species,val);
}

//  constrain species amount to a specified upper bound, units= moles,kg,m3
void GEMSEngine::set_species_upper_bound(const std::string& species, double val, const std::string& units)
{
    auto species_idx =gem.indexSpecies(species);
    if( units == "kg") {
        val/=species_molar_mass[species];
    }
    if( units == "m3") {
        val/=species_molar_volumes[species];
    }
    gem.setSpeciesUpperLimit(species,val);
}

// supresses a phase in GEM calculation
void GEMSEngine::supress_phase(const std::string& phase_name)
{
    supress_multiple_species(species_in_phase[phase_name]);
}

// supresses multiple phase in calculation as given in phase names list
void GEMSEngine::supress_multiple_phases(const std::vector<std::string>& phase_name_list)
{
    for( const auto& phase: phase_name_list) {
        supress_phase(phase);
    }
}

// supresses species in calculation
void GEMSEngine::supress_species(const std::string& species_name)
{
    set_species_lower_bound(species_name, 0);
    set_species_upper_bound(species_name, 1e-15);
}

// supresses multiple species in GEM calculation as given in species name list
void GEMSEngine::supress_multiple_species(const std::vector<std::string>& species_list)
{
    for( const auto& species: species_list) {
        supress_species(species);
    }
}

// activate supressed phase
void GEMSEngine::activate_phase(const std::string& phase_name)
{
    activate_multiple_species(species_in_phase[phase_name]);
}

// activate multiple supressed phases given in list
void GEMSEngine::activate_multiple_phases(const std::vector<std::string>& phase_name_list)
{
    for( const auto& phase: phase_name_list) {
        activate_phase(phase);
    }
}

// activate multiple supressed species given in the list
void GEMSEngine::activate_multiple_species(const std::vector<std::string>& species_list)
{
    for( const auto& species: species_list) {
        activate_species(species);
    }
}

// activate a supressed species in phase
void GEMSEngine::activate_species(const std::string& species_name)
{
    set_species_lower_bound(species_name, 0);
    set_species_upper_bound(species_name, 1e6);
}


// returns pH of the solution
double GEMSEngine::pH()
{
    return gem.pH();
}

// returns pE of the solution
double GEMSEngine::pE()
{
    return gem.pe();
}

// returns ionic strength of the solution
double GEMSEngine::ionic_strength()
{
    return gem.ionicStrength();
}

// returns volume of the system in m3
double GEMSEngine::system_volume()
{
    return gem.systemVolume();
}

// returns mass of the system in kg
double GEMSEngine::system_mass()
{
    return gem.systemMass();
}

// returns molar volume of phases in m3/mol
ValuesMap GEMSEngine::phases_molar_volume()
{
    ValuesMap phase_mvol;
    for(Index i = 0; i < nphases(); ++i) {
        phase_mvol[phase_names[i]] = gem.phaseMolarVolume(i);
        }
    return phase_mvol;
}

// returns saturation indices of phases
ValuesMap GEMSEngine::phase_sat_indices()
{
    return to_map(phase_names, gem.phaseSatIndices());
}

}
