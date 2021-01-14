// xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2018 Allan Leal, Dmitrii Kulik
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "ChemicalEngine.hpp"

// C++ includes
#include <chrono>
#include <iomanip>
#include <memory>

// GEMS3K includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#define NODEARRAYLEVEL
#include <GEMS3K/node.h>
#include <GEMS3K/gems3k_impex.h>

#define NODE  pimpl->node
#define CSD   pimpl->node->pCSD()

namespace xGEMS {

struct ChemicalEngine::Impl
{
    /// The TNode instance from ChemicalEngine
    std::unique_ptr<TNode> node;

    /// The options for ChemicalEngine
    ChemicalEngineOptions options;

    /// Molar masses of elements
    Vector elemMolMasses;

    /// Molar masses of species
    Vector molMasses;
    
    /// The formula matrix of the species
    Matrix formula_matrix;

    /// The elapsed time of the equilibrate method (in units of s)
    double elapsed_time = 0;

    /// Mole fractions of species in their respective phases 
    Vector molFractions;

    /// Mole fractions of species in their respective phases 
    Vector molalities;

    /// Activity coefficients of the species (in natural log scale)
    Vector ln_activity_coefficients;

    /// Activities of the species (in natural log scale)
    Vector ln_activities;

    /// Concentrations of the species (in natural log scale)
    Vector ln_concentrations;

    /// Chemical potentials of the species (in J/mol)
    Vector chemPotentials;

    /// Standard thermodynamic properties of species at T,P of the node
    // Vector stMolarGibbsEnergies;            // in J/mol
    // Vector stMolarEnthalpies;               // in J/mol
    // Vector stMolarEntropies;                // in J/K/mol
    // Vector stMolarHeatCapacitiesConstP;     // in J/K/mol
    // Vector stMolarVolumes;                  // in m3/mol

    /// molar properties of the phases .
    // Vector phMolarVolumes;             // (in units of m3/mol)
    // Vector phMolarGibbsEnergies;       // (in units of J/mol)
    // Vector phMolarEnthalpies;          // (in units of J/mol)
    // Vector phMolarEntropies;             // J/K/mol
    Vector phMolarHeatCapacitiesConstP;  // J/K/mol

    /// Densities of the phases (in units of kg/m3).
    Vector phDensities;

    /// Volumes of the phases (in units of m3).
    Vector phVolumes;
    
    /// Enthalpies of the phases (in units of J).
    Vector phEnthalpies;

     /// Entropies of the phases (in units of J/K).
    Vector phEntropies;

    /// Entropies of the phases (in units of J/K).
    Vector phHeatCapacitiesConstP;

    /// Mole amounts of the phases (in units of mol).
    Vector phAmounts;

    /// Masses of the phases (in units of kg).
    Vector phMasses;

    /// Saturation (stability) indexes of phases (in log10 scale).
    Vector phSatIndices;

    /// Json string with I/O node (DBR) object in JSON format
    std::string dbr_json;

    /// Write IPM, DCH and DBR files (in binary, txt or json mode)
    GEMS3KGenerator::IOModes io_mode = GEMS3KGenerator::f_key_value;


    /// Construct a default Impl instance
    Impl()
    {}
};

ChemicalEngine::ChemicalEngine()
: pimpl(new Impl())
{
}

ChemicalEngine::ChemicalEngine(std::string filename)
: pimpl(new Impl())
{
    initialize(filename);
}

ChemicalEngine::~ChemicalEngine()
{}

auto ChemicalEngine::initialize(std::string filename) -> void
{
    // Initialize the GEMS input files mode
    GEMS3KGenerator generator_data(filename.c_str());
    pimpl->io_mode =   generator_data.files_mode();

    // Allocate memory for the GEMS `node` member
    pimpl->node = std::unique_ptr<TNode>(new TNode);

    // Initialize the GEMS `node` member
    long int res = pimpl->node->GEM_init(filename.c_str());

    // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not initialize the ChemicalEngine object.\n"
                "There was a problem during memory allocation.");

    // Check if there was a file read error during node initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not initialize the ChemicalEngine object.\n"
                "Make sure the provided file exists relative to the working directory.");

    // Assemble the formula matrix of the species
    const Index E = numElements();
    const Index N = numSpecies();
    const Index P = numPhases();
    pimpl->formula_matrix = Matrix::Zero(E, N);
    for(Index i = 0; i < N; ++i)
        for(Index j = 0; j < E; ++j)
            pimpl->formula_matrix(j, i) = pimpl->node->DCaJI(i, j);

    // Allocate memory for vector members
    pimpl->elemMolMasses.resize(E);
    pimpl->molMasses.resize(N);
    pimpl->ln_activity_coefficients.resize(N);
    pimpl->ln_activities.resize(N);
    pimpl->ln_concentrations.resize(N);
    pimpl->chemPotentials.resize(N);
    pimpl->molFractions.resize(N);
    pimpl->molalities.resize(N);
    // pimpl->phMolarVolumes.resize(P);
    // pimpl->phMolarEnthalpies.resize(P);
    // pimpl->phMolarEntropies.resize(P);
    pimpl->phMolarHeatCapacitiesConstP.resize(P);
    pimpl->phDensities.resize(P);
    pimpl->phVolumes.resize(P);
    pimpl->phEnthalpies.resize(P);
    pimpl->phEntropies.resize(P);
    pimpl->phHeatCapacitiesConstP.resize(P);
    pimpl->phAmounts.resize(P);
    pimpl->phMasses.resize(P);
    pimpl->phSatIndices.resize(P);
    // pimpl->stMolarGibbsEnergies.resize(N);            // in J/mol
    // pimpl->stMolarEnthalpies.resize(N);               // in J/mol
    // pimpl->stMolarEntropies.resize(N);                // in J/K/mol
    // pimpl->stMolarHeatCapacitiesConstP.resize(N);     // in J/K/mol
    // pimpl->stMolarVolumes.resize(N);                  // in m3/mol
}


auto ChemicalEngine::initializeJstr(std::string& dch_json, std::string& ipm_json, std::string& dbr_json) -> void
{
    // pimpl->io_mode = "json";
    // Allocate memory for the GEMS `node` member
    pimpl->node = std::unique_ptr<TNode>(new TNode);

    // Initialize the GEMS `node` member
    long int res = pimpl->node->GEM_init( dch_json, ipm_json, dbr_json );

    // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not initialize the ChemicalEngine object from JSON strings. \n"
                " There was a problem during memory allocation.");

    // Check if there was a JSON parsing error during node initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not initialize the ChemicalEngine object from JSON strings.\n"
                "Make sure that these JSON strings are not empty and are in correct format.");

    // Assemble the formula matrix of the species
    const Index E = numElements();
    const Index N = numSpecies();
    const Index P = numPhases();
    pimpl->formula_matrix = Matrix::Zero(E, N);
    for(Index i = 0; i < N; ++i)
        for(Index j = 0; j < E; ++j)
            pimpl->formula_matrix(j, i) = pimpl->node->DCaJI(i, j);

    // Allocate memory for vector members
    pimpl->elemMolMasses.resize(E);
    pimpl->molMasses.resize(N);
    pimpl->ln_activity_coefficients.resize(N);
    pimpl->ln_activities.resize(N);
    pimpl->ln_concentrations.resize(N);
    pimpl->chemPotentials.resize(N);
    pimpl->molFractions.resize(N);
    pimpl->molalities.resize(N);
    // pimpl->phMolarVolumes.resize(P);
    // pimpl->phMolarEnthalpies.resize(P);
    // pimpl->phMolarEntropies.resize(P);
    pimpl->phMolarHeatCapacitiesConstP.resize(P);
    pimpl->phDensities.resize(P);
    pimpl->phVolumes.resize(P);
    pimpl->phEnthalpies.resize(P);
    pimpl->phEntropies.resize(P);
    pimpl->phHeatCapacitiesConstP.resize(P);
    pimpl->phAmounts.resize(P);
    pimpl->phMasses.resize(P);
    pimpl->phSatIndices.resize(P);
    // pimpl->stMolarGibbsEnergies.resize(N);            // in J/mol
    // pimpl->stMolarEnthalpies.resize(N);               // in J/mol
    // pimpl->stMolarEntropies.resize(N);                // in J/K/mol
    // pimpl->stMolarHeatCapacitiesConstP.resize(N);     // in J/K/mol
    // pimpl->stMolarVolumes.resize(N);                  // in m3/mol
}


auto ChemicalEngine::readDbrFile(std::string filename) -> void
{
    // Reads another dbr file with input system composition
    long int res = pimpl->node->GEM_read_dbr(filename.c_str(), pimpl->io_mode);

        // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not find the provided dbr file.\n"
                "There was a problem during memory allocation.");

    // Check if there was a file read error during node initialization
    if(res >= 1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not find the provided dbr file.\n"
                "Make sure the provided file path exists relative to the working directory.");
}

auto ChemicalEngine::readDbrJstr(std::string& dbr_json) -> void
{
    // Reads another dbr file with input system composition
    long int res = pimpl->node->GEM_read_dbr( dbr_json, true );
    // res = 1;    // Temporary plug

        // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not process the provided DBR JSON string.\n"
                "There was a problem during memory allocation.");

    // Check if there was a JSON string read error during node initialization
    if(res >= 1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not process the provided DBR JSON string.\n"
                "Make sure that this JSON string is not empty and is in correct format and parameter position.");
}

auto ChemicalEngine::writeDbrFile(std::string filename) -> void
{
    // Write current node into a dbr file into path filename
    pimpl->node->GEM_write_dbr(filename.c_str(), pimpl->io_mode);
}

auto ChemicalEngine::writeDbrJstr(std::string& dbr_json) -> void
{
    // Writes another node dbr with input system composition into json string dbr_json
    long int res = pimpl->node->GEM_write_dbr( dbr_json );

    // Check if there was a json write error
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not process the provided DBR into JSON string.\n"
                "Check the code.");
}

auto ChemicalEngine::numElements() const -> Index
{
    return pimpl->node->pCSD()->nIC;
}

auto ChemicalEngine::numSpecies() const -> Index
{
    return pimpl->node->pCSD()->nDC;
}

auto ChemicalEngine::numPhases() const -> Index
{
    return pimpl->node->pCSD()->nPH;
}

auto ChemicalEngine::numSpeciesInPhase(Index iphase) const -> Index
{
    return pimpl->node->pCSD()->nDCinPH[iphase];
}

auto ChemicalEngine::elementName(Index ielement) const -> std::string
{
    std::string name = pimpl->node->pCSD()->ICNL[ielement];
    if( name.length() > MaxICnameLength )
        name.resize( MaxICnameLength );
    return name;
}

auto ChemicalEngine::speciesName(Index ispecies) const -> std::string
{
    std::string name = pimpl->node->pCSD()->DCNL[ispecies];
       if( name.length() > MaxDCnameLength )
        name.resize( MaxDCnameLength );
    return name;
}

auto ChemicalEngine::speciesCharge(Index ispecies) const -> double
{
    MatrixConstRef W = formulaMatrix();
    auto chargeRow = numElements();
    return W(chargeRow -1 , ispecies);
}

auto ChemicalEngine::phaseName(Index iphase) const -> std::string
{
    std::string name = pimpl->node->pCSD()->PHNL[iphase];
       if( name.length() > MaxPHnameLength )
        name.resize( MaxPHnameLength );
    return name;
}

// These methods may be ambiguous, as in GEMS3K, species with the same name 
//   may occur in more than one condensed phase!
// An overload including phase name or index is needed! 
//
auto ChemicalEngine::setSpeciesUpperLimit(std::string name, double amount) -> void
{
    auto ispecies = indexSpecies(name);
    double bound = (amount < 0.0? 1e6: amount);
    pimpl->node->pCNode()->dul[ispecies] = bound;
}

auto ChemicalEngine::setSpeciesLowerLimit(std::string name, double amount) -> void
{
    auto ispecies = indexSpecies(name);
    double bound = (amount < 0.0? 0.0: amount);
    pimpl->node->pCNode()->dll[ispecies] = bound;
}

// Overload!
auto ChemicalEngine::setSpeciesUpperLimit(Index ispecies, double amount) -> void
{   
    double bound = (amount < 0.0? 1e6: amount);
    pimpl->node->pCNode()->dul[ispecies] = bound;
}

// Overload!
auto ChemicalEngine::setSpeciesLowerLimit(Index ispecies, double amount) -> void
{
    double bound = (amount < 0.0? 0.0: amount);
    pimpl->node->pCNode()->dll[ispecies] = bound;
}

auto ChemicalEngine::indexElement(std::string element) const -> Index
{
    const Index size = numElements();
    for(Index i = 0; i < size; ++i)
        if(elementName(i) == element)
            return i;
    return size; // in case that element name was not found 
}

auto ChemicalEngine::indexSpecies(std::string species) const -> Index
{
    const Index size = numSpecies();
    for(Index i = 0; i < size; ++i)
        if(speciesName(i) == species)
            return i;
    return size; // in case that species name was not found 
}

auto ChemicalEngine::indexSpeciesAll(std::string species) const -> VectorXi
{
    std::vector<Index> occur;
    occur.resize(0);
    const Index size = numSpecies();
    for(Index i = 0; i < size; ++i)
        if(speciesName(i) == species)
            occur.push_back(i);
    Index n_species_found = occur.size();
    VectorXi result;
    result.resize(0);
    if(n_species_found > 0)        
    {
        result.resize(n_species_found);
        for(Index is=0; is<n_species_found; is++)
            result(is) = occur.at(is);
    }
    return result;
}

// Index of phase searched by name
auto ChemicalEngine::indexPhase(std::string phase) const -> Index
{
    const Index size = numPhases();
    for(Index i = 0; i < size; ++i)
        if(phaseName(i) == phase)
            return i;
    return size; // in case that phase name was not found 
}

auto ChemicalEngine::indexPhaseAll(std::string phase) const -> VectorXi
{
    std::vector<Index> occur;
    occur.resize(0);
    const Index size = numPhases();
    for(Index i = 0; i < size; ++i)
        if(phaseName(i) == phase)
            occur.push_back(i);
    Index n_phases_found = occur.size();
    VectorXi result;
    result.resize(0);
    if(n_phases_found > 0)        
    {
        result.resize(n_phases_found);
        for(Index is=0; is<n_phases_found; is++)
            result(is) = occur.at(is);
    }
    return result;
}

auto ChemicalEngine::indexPhaseWithSpecies(Index ispecies) const -> Index
{
    const Index size = numPhases();
    Index counter = 0;
    for(Index i = 0; i < size; ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > ispecies) return i;
    }
    return size; // in case that species name was not found 
}

auto ChemicalEngine::indexFirstSpeciesInPhase(Index iphase) const -> Index
{
    Index counter = 0;
    for(Index i = 0; i < iphase; ++i)
        counter += numSpeciesInPhase(i);
    return counter;
}

auto ChemicalEngine::elementMolarMasses() const -> VectorConstRef
{
    for(Index i = 0; i < numElements(); ++i)
        pimpl->elemMolMasses[i] = pimpl->node->pCSD()->ICmm[i];
    return pimpl->elemMolMasses;
}

auto ChemicalEngine::speciesMolarMasses() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->molMasses[i] = pimpl->node->pCSD()->DCmm[i];
    return pimpl->molMasses;
}

auto ChemicalEngine::formulaMatrix() const -> MatrixConstRef
{
    return pimpl->formula_matrix;
}

auto ChemicalEngine::setOptions(const ChemicalEngineOptions& options) -> void
{
    pimpl->options = options;
}

auto ChemicalEngine::setWarmStart() -> void
{
    pimpl->options.warmstart = true;
}

auto ChemicalEngine::setColdStart() -> void
{
    pimpl->options.warmstart = false;
}

auto ChemicalEngine::setPT( double P, double T) const -> bool
{
    // Check if P, T are feasible
    if(pimpl->node->check_TP( T, P )==false)
        return true;
    // Set temperature and pressure
    pimpl->node->Set_TK(T);
    pimpl->node->Set_P(P);
    return false;    
}  

auto ChemicalEngine::setB( VectorConstRef b ) -> void
{
    // Set the mole amounts of elements
    for(Index i = 0; i < numElements(); ++i)
        pimpl->node->Set_bIC(i, b[i]);
}  

auto ChemicalEngine::reequilibrate( bool warmstart ) -> int
{
    pimpl->options.warmstart = warmstart;
    // Begin timing
    auto begin = std::chrono::high_resolution_clock::now();
    // Solve the equilibrium problem with gems
    pimpl->node->pCNode()->NodeStatusCH = 
       pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    auto valueOutputGem = pimpl->node->GEM_run(false);
    // Finish timing
    auto end = std::chrono::high_resolution_clock::now();
    // Set the elapsed time member
    pimpl->elapsed_time = std::chrono::duration<double>(end - begin).count();
    return valueOutputGem;
}

auto ChemicalEngine::equilibrate(double T, double P, VectorConstRef b) -> int
{
    // Begin timing
    auto begin = std::chrono::high_resolution_clock::now();

    // Set temperature and pressure
    pimpl->node->Set_TK(T);
    pimpl->node->Set_P(P);

    // Set the mole amounts of the elements
    for(Index i = 0; i < numElements(); ++i)
        pimpl->node->Set_bIC(i, b[i]);

    // Solve the equilibrium problem with gems
    pimpl->node->pCNode()->NodeStatusCH = 
       pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    auto valueOutputGem = pimpl->node->GEM_run(false);

    // Finish timing
    auto end = std::chrono::high_resolution_clock::now();

    // Set the elapsed time member
    pimpl->elapsed_time = std::chrono::duration<double>(end - begin).count();
    return valueOutputGem;
}

auto ChemicalEngine::setSpeciesAmounts(VectorConstRef n) -> void
{

    MatrixConstRef W = formulaMatrix();
    // Updates all the species amount
    long int jj;
    for(jj = 0; jj <numSpecies(); jj++)
        pimpl->node->pCNode()->xDC[jj] = n[jj];
    // Correction of the bIC vector
    long int ii;

    for( ii=0; ii<numElements(); ii++ )
        pimpl->node->pCNode()->bIC[ii] = 0.0;
    for(jj = 0; jj <numSpecies(); jj++)
        for( ii=0; ii<numElements(); ii++ )
            pimpl->node->pCNode()->bIC[ii] += speciesAmounts()[jj] * W(ii , jj);
}

// Set the species of name (string) with the new amount
// Caution: this may be ambiguous as in GEMS3K, species with the same name 
//   may occur in more than one condensed phase!
// An overload including phase name is needed! 
auto ChemicalEngine::setSpeciesAmount(std::string name, double amount) -> void
{

    MatrixConstRef W = formulaMatrix();
    auto ispecies = indexSpecies(name);
    // Updates the species amount
    pimpl->node->pCNode()->xDC[ispecies] = amount;
    // Correction of the bIC vector
    long int ii;
    long int jj;
    for( ii=0; ii<numElements(); ii++ )
        pimpl->node->pCNode()->bIC[ii] = 0.0;
    for(jj = 0; jj <numSpecies(); jj++)
        for( ii=0; ii<numElements(); ii++ )
            pimpl->node->pCNode()->bIC[ii] += speciesAmounts()[jj] * W(ii , jj);
}


auto ChemicalEngine::setSpeciesAmount(Index ispecies, double amount) -> void
{

    MatrixConstRef W = formulaMatrix();
    // Updates the species amount
    pimpl->node->pCNode()->xDC[ispecies] = amount;
    // Correction of the bIC vector
    long int ii;
    long int jj;
    for( ii=0; ii<numElements(); ii++ )
        pimpl->node->pCNode()->bIC[ii] = 0.0;
    for(jj = 0; jj <numSpecies(); jj++)
        for( ii=0; ii<numElements(); ii++ )
            pimpl->node->pCNode()->bIC[ii] += speciesAmounts()[jj] * W(ii , jj);
}

auto ChemicalEngine::converged() const -> bool
{
    const auto status = pimpl->node->pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA || status == OK_GEM_SIA;
}

auto ChemicalEngine::numIterations() const -> Index
{
    return pimpl->node->pCNode()->IterDone;
}

auto ChemicalEngine::elapsedTime() const -> double
{
    return pimpl->elapsed_time;
}

auto ChemicalEngine::temperature() const -> double
{
    return pimpl->node->Get_TK();
}

auto ChemicalEngine::pressure() const -> double
{
    return pimpl->node->Get_P();
}

auto ChemicalEngine::elementAmounts() const -> VectorConstRef
{
    return Vector::Map(pimpl->node->pCNode()->bIC, numElements());
}

auto ChemicalEngine::elementAmountsInPhase(Index iphase) const -> Vector
{
    MatrixConstRef W = formulaMatrix();
    VectorConstRef n = speciesAmounts();
    const Index first = indexFirstSpeciesInPhase(iphase);
    const Index size = numSpeciesInPhase(iphase);
    MatrixConstRef Wp = W.middleCols(first, size);
    VectorConstRef np = n.segment(first, size);
    Vector res = Wp * np;
    return res;
}

auto ChemicalEngine::elementAmountsInSpecies(VectorXiConstRef ispecies) const -> Vector
{
    MatrixConstRef W = formulaMatrix();
    VectorConstRef n = speciesAmounts();
    MatrixConstRef Wp = W(all, ispecies);
    VectorConstRef np = n(ispecies);
    Vector res = Wp * np;
    return res;
}

auto ChemicalEngine::speciesAmount(Index ispecies) const -> double
{
    return pimpl->node->pCNode()->xDC[ispecies];
}

auto ChemicalEngine::speciesAmount(std::string name) const -> double
{
    return pimpl->node->pCNode()->xDC[indexSpecies(name)];
}

auto ChemicalEngine::speciesAmounts() const -> VectorConstRef
{
    return Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
}

// Aquatic systems only (assuming aqueous phase is the first one and H2O-solvent 
// is the last species in it)
auto ChemicalEngine::speciesMolalities() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->molalities[i] = pimpl->node->Get_cDC(i);
    if(numSpecies() > numSpeciesInPhase(0)) 
    {
        Vector amountSp = Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
        Index H2Oindex = numSpeciesInPhase(0)-1;
        double H2Omass = pimpl->node->Get_nDC(H2Oindex)*pimpl->node->DCmm(H2Oindex); // in kg
        if( H2Omass != 0.0 )
            for(Index i = H2Oindex+1; i < numSpecies(); ++i)
                pimpl->molalities[i] = speciesAmounts()(i)/H2Omass; 
    }               
    return pimpl->molalities;
}

auto ChemicalEngine::moleFractions() const -> VectorConstRef
{
    Vector amountSp = Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
    Index counter = 0;
    for(Index k = 0; k < numPhases(); ++k)
    {
        Index numSpInPh = numSpeciesInPhase(k);
        double amountPh = phaseAmounts()(k);
        if(numSpInPh == 1 )
        {
            if(amountPh > 0.0)
                pimpl->molFractions(counter) = 1.0;
            else 
                pimpl->molFractions(counter) = 0.0; 
            counter++;   
            continue;
        }
        for(Index i = 0; i < numSpInPh; ++i)
        {
            if(amountPh != 0.0)
                pimpl->molFractions(counter) = amountSp(counter)/amountPh;
            else
                pimpl->molFractions(counter) = 0.0;
            counter++;    
        }
    }
    return pimpl->molFractions;
}

auto ChemicalEngine::lnActivityCoefficients() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->ln_activity_coefficients[i] = std::log(pimpl->node->Get_gDC(i));
    return pimpl->ln_activity_coefficients;
}

auto ChemicalEngine::lnActivities() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->ln_activities[i] = std::log(pimpl->node->Get_aDC(i, true));
    return pimpl->ln_activities;
}

auto ChemicalEngine::lnConcentrations() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->ln_concentrations[i] = std::log(pimpl->node->Get_cDC(i));
    return pimpl->ln_concentrations;
}

auto ChemicalEngine::chemicalPotentials() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->chemPotentials[i] = pimpl->node->Get_muDC(i, false);
    return pimpl->chemPotentials;
    
    return Vector{};
}

// Implemented on Oct 1, 2020
// auto ChemicalEngine::standardMolarGibbsEnergies() const -> VectorConstRef
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     for(Index i = 0; i < numSpecies(); ++i)
//     {
//         auto idc = pimpl->node->DC_xDB_to_xCH( i );
//         pimpl->stMolarGibbsEnergies[i] = pimpl->node->DC_G0(idc, P, TK,  false);
//     }    
//     return pimpl->stMolarGibbsEnergies;
// }

// Uses curent node TK and P
auto ChemicalEngine::standardMolarGibbsEnergy(Index ispecies) const -> double
{
    double TK = pimpl->node->cTK();
    double P = pimpl->node->cP();
    auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
    auto stMolarG = pimpl->node->DC_G0(idc, P, TK,  false);  
    return stMolarG;
}

// Implemented on Oct 1, 2020
auto ChemicalEngine::standardMolarEnthalpy(Index ispecies) const -> double
{
    double TK = pimpl->node->cTK();
    double P = pimpl->node->cP();
    auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
    auto stMolarH = pimpl->node->DC_H0(idc, P, TK);   
    return stMolarH;
}

// Implemented on Oct 1, 2020
auto ChemicalEngine::standardMolarVolume(Index ispecies) const -> double
{
    double TK = pimpl->node->cTK();
    double P = pimpl->node->cP();
    auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
    auto stMolarV = pimpl->node->DC_V0(idc, P, TK);   
    return stMolarV;
}

// Implemented on Oct 1, 2020
auto ChemicalEngine::standardMolarEntropy(Index ispecies) const -> double
{
    double TK = pimpl->node->cTK();
    double P = pimpl->node->cP();
    auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
    auto stMolarS = pimpl->node->DC_S0(idc, P, TK); 
    return stMolarS;
}

// TBD
// auto ChemicalEngine::standardMolarInternalEnergy(Index ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarU = pimpl->node->DC_U0(idc, P, TK); 
//     return stMolarU;
// }

// TBD
// auto ChemicalEngine::standardMolarHelmholtzEnergy(Index ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarA = pimpl->node->DC_A0(idc, P, TK); 
//     return stMolarA;
// }

// Implemented on Oct 1, 2020
auto ChemicalEngine::standardMolarHeatCapacityConstP(Index ispecies) const -> double
{
    double TK = pimpl->node->cTK();
    double P = pimpl->node->cP();
    auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
    auto stMolarCp = pimpl->node->DC_Cp0(idc, P, TK); 
    return stMolarCp;
}

// TBD
// auto ChemicalEngine::standardMolarHeatCapacityConstV(Index ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarCv = pimpl->node->DC_Cp0(idc, P, TK) - R_CONSTANT; 
//     return stMolarCv;
// }

auto ChemicalEngine::phaseMolarGibbsEnergy(Index iphase) const -> double
{
    double phMolarGibbsEnergy = 0.0;
    if( pimpl->node->Ph_Mole(iphase) > 1e-15 )
        phMolarGibbsEnergy = pimpl->node->Ph_GibbsEnergy(iphase)/pimpl->node->Ph_Mole(iphase);
    return phMolarGibbsEnergy;
}

auto ChemicalEngine::phaseMolarEnthalpy(Index iphase) const -> double
{
    double phMolarEnthalpy = 0.0;
    if( pimpl->node->Ph_Mole(iphase) > 1e-15 )
        phMolarEnthalpy = pimpl->node->Ph_Enthalpy(iphase)/pimpl->node->Ph_Mole(iphase);
    return phMolarEnthalpy;
}

auto ChemicalEngine::phaseMolarVolume(Index iphase) const -> double
{
     double phMolarVolume = 0.0;
    if( pimpl->node->Ph_Mole(iphase) > 1e-15 )
        phMolarVolume = pimpl->node->Ph_Volume(iphase)/pimpl->node->Ph_Mole(iphase);
    return phMolarVolume;
}

auto ChemicalEngine::phaseMolarEntropy(Index iphase) const -> double
{
    double phMolarEntropy = 0.0;
    if( pimpl->node->Ph_Mole(iphase) > 1e-15 )
        phMolarEntropy = pimpl->node->Ph_Entropy(iphase)/pimpl->node->Ph_Mole(iphase);
    return phMolarEntropy;
}

// TBD
// auto ChemicalEngine::phaseMolarInternalEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

// TBD
// auto ChemicalEngine::phaseMolarHelmholtzEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

auto ChemicalEngine::phaseMolarHeatCapacityConstP(Index iphase) const -> double
{
     double phMolarCp = 0.0;
    if( pimpl->node->Ph_Mole(iphase) > 1e-15 )
        phMolarCp = pimpl->node->Ph_HeatCapacityCp(iphase)/pimpl->node->Ph_Mole(iphase);
    return phMolarCp;
}

// TBD
// auto ChemicalEngine::phaseMolarHeatCapacityConstV(Index iphase) const -> double
// {
//     return 0.0;
// }

auto ChemicalEngine::phaseSpecificGibbsEnergy(Index iphase) const -> double
{
     double phSpecGibbsEnergy = 0.0;
    if( pimpl->node->Ph_Mass(iphase) > 1e-15 )
        phSpecGibbsEnergy = pimpl->node->Ph_GibbsEnergy(iphase)/pimpl->node->Ph_Mass(iphase);
    return phSpecGibbsEnergy;
}

auto ChemicalEngine::phaseSpecificEnthalpy(Index iphase) const -> double
{
    double phSpecEnthalpy = 0.0;
    if( pimpl->node->Ph_Mass(iphase) > 1e-15 )
        phSpecEnthalpy = pimpl->node->Ph_Enthalpy(iphase)/pimpl->node->Ph_Mass(iphase);
    return phSpecEnthalpy;
}

auto ChemicalEngine::phaseSpecificVolume(Index iphase) const -> double
{
     double phSpecVolume = 0.0;
    if( pimpl->node->Ph_Mass(iphase) > 1e-15 )
        phSpecVolume = pimpl->node->Ph_Volume(iphase)/pimpl->node->Ph_Mass(iphase);
    return phSpecVolume;
}

auto ChemicalEngine::phaseSpecificEntropy(Index iphase) const -> double
{
    double phSpecEntropy = 0.0;
    if( pimpl->node->Ph_Mass(iphase) > 1e-15 )
        phSpecEntropy = pimpl->node->Ph_Entropy(iphase)/pimpl->node->Ph_Mass(iphase);
    return phSpecEntropy;
}

// TBD
// auto ChemicalEngine::phaseSpecificInternalEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

// TBD
// auto ChemicalEngine::phaseSpecificHelmholtzEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

auto ChemicalEngine::phaseSpecificHeatCapacityConstP(Index iphase) const -> double
{
    double phSpecCp = 0.0;
    if( pimpl->node->Ph_Mass(iphase) > 1e-15 )
        phSpecCp = pimpl->node->Ph_HeatCapacityCp(iphase)/pimpl->node->Ph_Mass(iphase);
    return phSpecCp;
}

// auto ChemicalEngine::phaseSpecificHeatCapacityConstV(Index iphase) const -> double
// {
//     return 0.0;
// }

auto ChemicalEngine::phaseDensities() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
    {    
        if(pimpl->node->Ph_Volume(i) != 0.0)
            pimpl->phDensities[i] = pimpl->node->Ph_Mass(i)/pimpl->node->Ph_Volume(i);
        else 
            pimpl->phDensities[i] = 0.0;
    }
    return pimpl->phDensities;
}

auto ChemicalEngine::phaseMasses() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phMasses[i] = pimpl->node->Ph_Mass(i);
    return pimpl->phMasses;
}

auto ChemicalEngine::phaseAmounts() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phAmounts[i] = pimpl->node->Ph_Moles(i);
    return pimpl->phAmounts;
}

auto ChemicalEngine::phaseVolumes() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phVolumes[i] = pimpl->node->Ph_Volume(i);
    return pimpl->phVolumes;
}

auto ChemicalEngine::phaseEnthalpies() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phEnthalpies[i] = pimpl->node->Ph_Enthalpy(i);
    return pimpl->phEnthalpies;
}

auto ChemicalEngine::phaseEntropies() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phEntropies[i] = pimpl->node->Ph_Entropy(i);
    return pimpl->phEntropies;
}

auto ChemicalEngine::phaseHeatCapacitiesConstP() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phHeatCapacitiesConstP[i] = pimpl->node->Ph_HeatCapacityCp(i);
    return pimpl->phHeatCapacitiesConstP;
}

auto ChemicalEngine::phaseDensity(Index iphase) const -> double
{
    if(pimpl->node->Ph_Volume(iphase) > 1e-15)
        return pimpl->node->Ph_Mass(iphase)/pimpl->node->Ph_Volume(iphase);
    else return 0.0;    
}

auto ChemicalEngine::phaseMass(Index iphase) const -> double
{
    return pimpl->node->Ph_Mass(iphase);
}

auto ChemicalEngine::phaseAmount(Index iphase) const -> double
{
    return pimpl->node->Ph_Moles(iphase);
}

auto ChemicalEngine::phaseVolume(Index iphase) const -> double
{
    return pimpl->node->Ph_Volume(iphase);
}

auto ChemicalEngine::phaseEnthalpy(Index iphase) const -> double
{
    return pimpl->node->Ph_Enthalpy(iphase);
}

auto ChemicalEngine::phaseEntropy(Index iphase) const -> double
{
    return pimpl->node->Ph_Entropy(iphase);
}

auto ChemicalEngine::phaseHeatCapacityConstP(Index iphase) const -> double
{
    return pimpl->node->Ph_HeatCapacityCp(iphase);
}

auto ChemicalEngine::phaseSatIndices() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phSatIndices[i] = pimpl->node->Ph_SatInd(i); 
    return pimpl->phSatIndices;
}

auto ChemicalEngine::phaseSatIndex(Index iphase) const -> double
{
    return pimpl->node->Ph_SatInd(iphase);
}

auto ChemicalEngine::systemVolume() const -> double
{
    return pimpl->node->pCNode()->Vs;
    // return pimpl->node->cVs();
}

auto ChemicalEngine::systemMass() const -> double
{
    //return pimpl->node->cMs();
    return pimpl->node->pCNode()->Ms;
}

auto ChemicalEngine::ionicStrength() const -> double
{
    return pimpl->node->pCNode()->IC;
}

auto ChemicalEngine::pH() const -> double
{
    return pimpl->node->pCNode()->pH;
}

auto ChemicalEngine::pe() const -> double
{
    return pimpl->node->pCNode()->pe;
}

auto ChemicalEngine::Eh() const -> double
{
    return pimpl->node->pCNode()->Eh;
}

auto ChemicalEngine::systemGibbsEnergy() const -> double
{
    return pimpl->node->pCNode()->Gs;   // *RT to get it in J
}

auto ChemicalEngine::systemEnthalpy() const -> double   // in J
{
    // enth = pimpl->node->phEnthalpies.sum(); // needs to be computed first anyway
    auto enth = 0.0;
    for(Index i = 0; i < numPhases(); ++i)
        enth += pimpl->node->Ph_Enthalpy(i);
    return enth;  // in J 
    // enth /= pimpl->node->pCNode()->Ms; in  kg
    // return pimpl->node->pCNode()->Hs;  // Not computed in GEMS3K so far!
}

auto ChemicalEngine::systemEntropy() const -> double
{
    return 0.0;
}

auto ChemicalEngine::systemHeatCapacityConstP() const -> double
{
    return 0.0;
}

auto operator<<(std::ostream& out, const ChemicalEngine& state) -> std::ostream&
{
    const double T = state.temperature();
    const double P = state.pressure();
    VectorConstRef n = state.speciesAmounts();
    const Vector activity_coeffs = state.lnActivityCoefficients().array().exp();
    const Vector activities = state.lnActivities().array() * 0.4343; // .exp();
    const Vector chemical_potentials = state.chemicalPotentials().array() / 1000.;
    const Vector concentrations = state.lnConcentrations().array().exp();
//    const Vector molalities = state.speciesMolalities().array();
    const Vector molfractions = state.moleFractions().array();
//    const Vector satindices = state.phaseSatIndices().array();

    const Index num_phases = state.numPhases();
    const Index bar_size = std::max(Index(9), num_phases + 2) * 25;
    const std::string bar1(bar_size, '=');
    const std::string bar2(bar_size, '-');

    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Temperature[K]";
    out << std::left << std::setw(25) << "Temperature[C]";
    out << std::left << std::setw(25) << "Pressure[MPa]";
    out << std::endl << bar2 << std::endl;

    out << std::left << std::setw(25) << T;
    out << std::left << std::setw(25) << T - 273.15;
    out << std::left << std::setw(25) << P * 1e-6;
    out << std::endl;

    // Set output in scientific notation
    auto flags = out.flags();
    out << std::setprecision(6);

    // Output the table of the element-related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Element";
    out << std::left << std::setw(25) << "Amount[mol]";
    for(Index i = 0; i < state.numPhases(); ++i)
        out << std::left << std::setw(25) << state.phaseName(i) + "[mol]";
    out << std::endl;
    out << bar2 << std::endl;
    for(Index i = 0; i < state.numElements(); ++i)
    {
        out << std::left << std::setw(25) << state.elementName(i);
        out << std::left << std::setw(25) << state.elementAmounts()[i];
        for(Index j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.elementAmountsInPhase(j)[i];
        out << std::endl;
    }
    out << bar2 << std::endl;
    out << std::left << std::setw(25) << "PhaseAmount[mol]";
    out << std::left << std::setw(25) << 0.0;
    for(Index j = 0; j < state.numPhases(); ++j)
        out << std::left << std::setw(25) << state.phaseAmounts()[j];
    out << std::endl;   
    out << std::left << std::setw(25) << "PhaseMass[kg]";
    out << std::left << std::setw(25) << 0.0;
    for(Index j = 0; j < state.numPhases(); ++j)
        out << std::left << std::setw(25) << state.phaseMasses()[j];
    out << std::endl;        
    out << std::left << std::setw(25) << "PhaseVolume[m^3]";
    out << std::left << std::setw(25) << 0.0;
    for(Index j = 0; j < state.numPhases(); ++j)
        out << std::left << std::setw(25) << state.phaseVolumes()[j];
    out << std::endl;  
    out << std::left << std::setw(25) << "PhaseDensity[kg/m^3]";
    out << std::left << std::setw(25) << 0.0;
    for(Index j = 0; j < state.numPhases(); ++j)
        out << std::left << std::setw(25) << state.phaseDensities()[j];
    out << std::endl;    
    out << std::left << std::setw(25) << "PhaseSatIndex[lg]";
    out << std::left << std::setw(25) << 0.0;
    for(Index j = 0; j < state.numPhases(); ++j)
//        out << std::left << std::setw(25) << satindices[j];    
        out << std::left << std::setw(25) << state.phaseSatIndices()[j];
    out << std::endl;
    out << std::endl;            

    // Output the table of the species-related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Species";
    out << std::left << std::setw(25) << "Amount[mol]";
    out << std::left << std::setw(25) << "Concentration[-]";
    out << std::left << std::setw(25) << "ActivityCoeff[-]";
 //   out << std::left << std::setw(25) << "Molality[mol/kgw]";
    out << std::left << std::setw(25) << "MoleFraction[-]";
    out << std::left << std::setw(25) << "ChemPotential[kJ/mol]";
    out << std::left << std::setw(25) << "Activity[lg]";
    out << std::endl;
    out << bar2 << std::endl;
    for(Index i = 0; i < state.numSpecies(); ++i)
    {
        out << std::left << std::setw(25) << state.speciesName(i);
        out << std::left << std::setw(25) << n[i];
        out << std::left << std::setw(25) << concentrations[i];
        out << std::left << std::setw(25) << activity_coeffs[i];
 //       out << std::left << std::setw(25) << molalities[i];
        out << std::left << std::setw(25) << molfractions[i];
        out << std::left << std::setw(25) << chemical_potentials[i];
        out << std::left << std::setw(25) << activities[i];
        out << std::endl;
    }

    // Output the table of the aqueous phase related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Mass[kg]";
    out << std::left << std::setw(25) << "Volume[m^3]";
    out << std::left << std::setw(25) << "G[norm]";
    out << std::left << std::setw(25) << "H[kJ]";
    out << std::left << std::setw(25) << "IonicStrength[molal]";
    out << std::left << std::setw(25) << "pH";
    out << std::left << std::setw(25) << "pE";
    out << std::left << std::setw(25) << "Eh[V]";
    out << std::endl << bar2 << std::endl;
    out << std::left << std::setw(25) << state.systemMass();
    out << std::left << std::setw(25) << state.systemVolume();
    out << std::left << std::setw(25) << state.systemGibbsEnergy();
    out << std::left << std::setw(25) << state.systemEnthalpy()/1000.0;
    out << std::left << std::setw(25) << state.ionicStrength();
    out << std::left << std::setw(25) << state.pH();
    out << std::left << std::setw(25) << state.pe();
    out << std::left << std::setw(25) << state.Eh();
    out << std::endl << bar1 << std::endl;

    // Recover the previous state of `out`
    out.flags(flags);

    return out;
}

} // namespace xGEMS
