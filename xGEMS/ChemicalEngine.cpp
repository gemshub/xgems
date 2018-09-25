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

    /// The activity coefficients of the species (in natural log scale)
    Vector ln_activity_coefficients;

    /// The activities of the species (in natural log scale)
    Vector ln_activities;

    /// The concentrations of the species (in natural log scale)
    Vector ln_concentrations;

    /// The chemical potentials of the species (in J/mol)
    Vector chemPotentials;

    /// molar volumes of the phases (in units of m3/mol).
    Vector phMolarVolumes;

    /// the densities of the phases (in units of kg/m3).
    Vector phDensities;

    /// volumes of the phases (in units of m3).
    Vector phVolumes;

    /// the molar amounts of the phases (in units of mol).
    Vector phAmounts;

    /// the masses of the phases (in units of kg).
    Vector phMasses;

    /// the saturation (stability) indexes of phases (in log10 scale).
    Vector phSatIndex;

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
    // Allocate memory for the GEMS `node` member
    pimpl->node = std::unique_ptr<TNode>(new TNode);

    // Initialize the GEMS `node` member
    const auto res = pimpl->node->GEM_init(filename.c_str());

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
    pimpl->phMolarVolumes.resize(P);
    pimpl->phDensities.resize(P);
    pimpl->phVolumes.resize(P);
    pimpl->phAmounts.resize(P);
    pimpl->phMasses.resize(P);
    pimpl->phSatIndex.resize(P);

}

auto ChemicalEngine::readDbrFile(std::string filename) -> void
{
    // Reads another dbr file with input system composition
    const auto res = pimpl->node->GEM_read_dbr(filename.c_str());

        // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not find the provided dbr file.\n"
                "There was a problem during memory allocation.");

    // Check if there was a file read error during node initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
            "Could not find the provided dbr file.\n"
                "Make sure the provided file exists relative to the working directory.");
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
    return name;
}

auto ChemicalEngine::speciesName(Index ispecies) const -> std::string
{
    std::string name = pimpl->node->pCSD()->DCNL[ispecies];
    return name;
}

auto ChemicalEngine::phaseName(Index iphase) const -> std::string
{
    std::string name = pimpl->node->pCSD()->PHNL[iphase];
    return name;
}

auto ChemicalEngine::indexElement(std::string element) const -> Index
{
    const Index size = numElements();
    for(Index i = 0; i < size; ++i)
        if(elementName(i) == element)
            return i;
    return size;
}

auto ChemicalEngine::indexSpecies(std::string species) const -> Index
{
    const Index size = numSpecies();
    for(Index i = 0; i < size; ++i)
        if(speciesName(i) == species)
            return i;
    return size;
}

auto ChemicalEngine::indexPhase(std::string phase) const -> Index
{
    const Index size = numPhases();
    for(Index i = 0; i < size; ++i)
        if(phaseName(i) == phase)
            return i;
    return size;
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
    return size;
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

auto ChemicalEngine::equilibrate(double T, double P, VectorConstRef b) -> void
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
    pimpl->node->GEM_run(false);

    // Finish timing
    auto end = std::chrono::high_resolution_clock::now();

    // Set the elapsed time member
    pimpl->elapsed_time = std::chrono::duration<double>(end - begin).count();
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

auto ChemicalEngine::speciesAmounts() const -> VectorConstRef
{
    return Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
}

auto ChemicalEngine::speciesMolalities() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->molalities[i] = pimpl->node->Get_aDC(i, true);
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
            if(amountPh > 0.0)
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

auto ChemicalEngine::standardPartialMolarGibbsEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarEnthalpies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarVolumes() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarEntropies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarInternalEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarHelmholtzEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarHeatCapacitiesConstP() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::standardPartialMolarHeatCapacitiesConstV() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarGibbsEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarEnthalpies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarVolumes() const -> VectorConstRef
{
/*
    VectorConstRef n = speciesAmounts();

    for(Index iphase = 0; iphase < numPhases(); ++iphase)
    {
        const Index first = indexFirstSpeciesInPhase(iphase);
        const Index size  = numSpeciesInPhase(iphase);
        VectorConstRef np = n.segment(first, size);  

        for (Index ispecies = 0; ispecies < size; ++ispecies)
        {
            // to m3/mol
            pimpl->phMolarVolumes[iphase] += np[ispecies] * pimpl->node->DC_V0(ispecies, temperature(), pressure()) * 1e-5;
        }
    }
    return pimpl->phMolarVolumes;
*/
    for(Index i = 0; i < numPhases(); ++i)
    {
        pimpl->phMolarVolumes[i] = 0.0;
        if( pimpl->node->Ph_Mole(i) )
            pimpl->phMolarVolumes[i] = pimpl->node->Ph_Volume(i)/pimpl->node->Ph_Mole(i);
    }
    return pimpl->phMolarVolumes;
}

auto ChemicalEngine::phaseMolarEntropies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarInternalEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarHelmholtzEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarHeatCapacitiesConstP() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseMolarHeatCapacitiesConstV() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificGibbsEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificEnthalpies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificVolumes() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificEntropies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificInternalEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificHelmholtzEnergies() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificHeatCapacitiesConstP() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseSpecificHeatCapacitiesConstV() const -> VectorConstRef
{
    return Vector{};
}

auto ChemicalEngine::phaseDensities() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phDensities[i] = pimpl->node->Ph_Mass(i)/pimpl->node->Ph_Volume(i);
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
        pimpl->phAmounts[i] = pimpl->node->Ph_Mole(i);
    return pimpl->phAmounts;
}

auto ChemicalEngine::phaseVolumes() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phVolumes[i] = pimpl->node->Ph_Volume(i);
    return pimpl->phVolumes;
}

auto ChemicalEngine::phaseSatIndices() const -> VectorConstRef
{
    for(Index i = 0; i < numPhases(); ++i)
        pimpl->phSatIndex[i] = pimpl->node->Ph_SatInd(i);
    return pimpl->phSatIndex;
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
    return pimpl->node->pCNode()->Gs;
}

auto ChemicalEngine::systemEnthalpy() const -> double
{
    return pimpl->node->pCNode()->Hs;
}

auto operator<<(std::ostream& out, const ChemicalEngine& state) -> std::ostream&
{
    const double T = state.temperature();
    const double P = state.pressure();
    VectorConstRef n = state.speciesAmounts();
    const Vector activity_coeffs = state.lnActivityCoefficients().array().exp();
    const Vector activities = state.lnActivities().array().exp();
    const Vector chemical_potentials = state.chemicalPotentials().array();
    const Vector concentrations = state.lnConcentrations().array().exp();
    const Vector molalities = state.speciesMolalities().array();

    const Index num_phases = state.numPhases();
    const Index bar_size = std::max(Index(9), num_phases + 2) * 25;
    const std::string bar1(bar_size, '=');
    const std::string bar2(bar_size, '-');

    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Temperature [K]";
    out << std::left << std::setw(25) << "Temperature [C]";
    out << std::left << std::setw(25) << "Pressure [MPa]";
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
    out << std::left << std::setw(25) << "Amount [mol]";
    for(Index i = 0; i < state.numPhases(); ++i)
        out << std::left << std::setw(25) << state.phaseName(i) + " [mol]";
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

    // Output the table of the species-related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Species";
    out << std::left << std::setw(25) << "Amount [mol]";
    out << std::left << std::setw(25) << "Activity Coeff. [-]";
    out << std::left << std::setw(25) << "Activity [-]";
    out << std::left << std::setw(25) << "Chem.Potential [kJ/mol]";
    out << std::left << std::setw(25) << "Concentration [-]";
    out << std::left << std::setw(25) << "Molality [mol/kgw]";
    out << std::endl;
    out << bar2 << std::endl;
    for(Index i = 0; i < state.numSpecies(); ++i)
    {
        out << std::left << std::setw(25) << state.speciesName(i);
        out << std::left << std::setw(25) << n[i];
        out << std::left << std::setw(25) << activity_coeffs[i];
        out << std::left << std::setw(25) << activities[i];
        out << std::left << std::setw(25) << chemical_potentials[i]/1000.0;
        out << std::left << std::setw(25) << concentrations[i];
        out << std::left << std::setw(25) << molalities[i];
        out << std::endl;
    }

    // Output the table of the aqueous phase related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "M [kg]";
    out << std::left << std::setw(25) << "V [m^3]";
    out << std::left << std::setw(25) << "G [norm]";
    out << std::left << std::setw(25) << "H [kJ]";
    out << std::left << std::setw(25) << "Ionic Strength [molal]";
    out << std::left << std::setw(25) << "pH";
    out << std::left << std::setw(25) << "pE";
    out << std::left << std::setw(25) << "Eh [V]";
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
