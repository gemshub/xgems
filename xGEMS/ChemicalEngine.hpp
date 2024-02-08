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

#pragma once

// C++ includes
#include <iostream>
#include <memory>
#include <string>

// xGEMS includes
#include <xGEMS/Index.hpp>
#include <xGEMS/Eigen.hpp>

namespace xGEMS {

/// Update loggers settings
/// @param use_cout:      show/hide logging to stdout
///        logfile_name:  add logging to rotating file name (hide if empty)
///        log_level:     set login level for all loggers
void update_loggers(bool use_cout, const std::string& logfile_name, size_t log_level);

/// A type that describes the options for ChemicalEngine
struct ChemicalEngineOptions
{
    /// The flag that indicates if smart start initial approximation is used
    bool warmstart = true;
};

/// A wrapper class for ChemicalEngine code
class ChemicalEngine
{
public:
    /// Construct a default ChemicalEngine object.
    ChemicalEngine();

    /// Construct a ChemicalEngine object from a GEM-Selektor project file.
    /// @param filename The name of the file containing the definition of the chemical system
    ChemicalEngine(std::string filename);

    /// Destroy this ChemicalEngine instance.
    virtual ~ChemicalEngine();

    /// Construct a copy of a ChemicalEngine object.
    ChemicalEngine(const ChemicalEngine& other) = delete;

    /// Assign another ChemicalEngine object to this.
    auto operator=(ChemicalEngine other) -> ChemicalEngine& = delete;

    /// Reallocates chemical engine vectors and matrix upon node initialization from files or json
    /// Used only inside of initialize() or initializeFromJsonStrings()
    auto reallocateEngineArrays() -> void;

    /// Initialize the ChemicalEngine object from a GEM-Selektor project file.
    /// @param filename The name of the file containing the definition of the chemical system
    auto initialize(std::string filename) -> void;
    
    /// Initialize the ChemicalEngine object from JSON strings for dch, ipm and dbr exported from GEM-Selektor \ .
    /// Normally to be used after ChemicalEngine();
    /// @param dch_json The json string containing the definition of the chemical system
    /// @param ipm_json The json string containing the parameters and settings for GEMS3K IPM-3 algorithm
    /// @param dbr_json The json string containing the input node composition of the chemical system
    auto initializeFromJsonStrings( std::string dch_json, std::string ipm_json, std::string dbr_json) -> void;

    /// Reads another DBR file (with input system composition, T,P etc.) \ . The DBR file must be compatible with
    /// the currently loaded IPM and DCH files (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFromFile(std::string filename) -> void;

    /// Reads another DBR object (with input system composition, T,P etc.) from JSON string \ . The DBR file 
    /// must be compatible with the currently loaded IPM and DCH objects (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFromJsonString( std::string dbr_json) -> void;

    /// Writes a DBR file (normally after some changes via API and GEM calculation). 
    /// @param Null-terminated (C) string containing a full path to the output DBR disk file.
    auto writeDbrToFile(std::string filename) -> void;

    /// Returns a DBR object (normally after some changes via API and GEM calculation) as JSON string. \ .
    /// In case of error raises an exception or returns empty string. 
    auto writeDbrToJsonString() -> const std::string;

    /// Return the number of elements.
    auto numElements() const -> Index;

    /// Return the number of species.
    auto numSpecies() const -> Index;

    /// Return the number of phases.
    auto numPhases() const -> Index;

    /// Return the number of species in a phase or 0, if the phase was not found
    /// @param iphase The index of the phase.
    auto numSpeciesInPhase(Index iphase) const -> Index;

    /// Return the name of an element.
    /// @param ielement The index of the element.
    /// In case of error raises an exception or returns empty string. 
    auto elementName(Index ielement) const -> std::string;

    /// Return the name of a species.
    /// @param ispecies The index of the species.
    /// In case of error raises an exception or returns empty string. 
    auto speciesName(Index ispecies) const -> std::string;

    /// Return charge of a species.
    auto speciesCharge(Index ispecies) const -> double;

    /// Return the name of a phase
    /// @param iphase The index of the phase.
    /// In case of error returns null string. 
    auto phaseName(Index iphase) const -> std::string;

    /// Return the index of an @param element name, or number of elements, if not found.
    /// @param Index The index of the element.
    auto indexElement(std::string element) const -> Index;

    /// Return the index of a @param species name, or number of species, if not found.
    /// @param Index The index of the species.
    auto indexSpecies(std::string species) const -> Index;

    /// Return the vector of all indices of a @param species name, number of species found is VectorXi.size().
    /// @param VectorXi contains all indices of the species found in the species name list in this system.
    auto indexSpeciesAll(std::string species) const -> VectorXi;

    /// Return the index of a @param phase name, or number of phases, if not found.
    /// @param Index The index of the phase.
    auto indexPhase(std::string phase) const -> Index;

    /// Return the vector of all indices of a @param phase name, number of indices found is VectorXi.size().
    /// @param VectorXi contains all indices of the phases found in the phase name list in this system.
    auto indexPhaseAll(std::string phase) const -> VectorXi;

    /// Return the index of the phase containing a @param species name (ambiguous).
    /// @param Index The index of the species.
    auto indexPhaseWithSpecies(Index ispecies) const -> Index;

    /// Return the @param Index of the first species in a phase.
    /// @param iphase The index of the phase.
    auto indexFirstSpeciesInPhase(Index iphase) const -> Index;

    /// Return the molar masses of the elements (in units of kg/mol)
    auto elementMolarMasses() const -> VectorConstRef;

    /// Return the molar masses of the species (in units of kg/mol)
    auto speciesMolarMasses() const -> VectorConstRef;

    /// Return the formula matrix of elements in species (in moles per mole species)
    /// (rows: elements, cols: species).
    auto formulaMatrix() const -> MatrixConstRef;

    /// Sets the whole speciation vector xDC to amounts provided in n and resets bIC bector accordingly
    ///  (often used for input initial guess of GEM solution in reactive transport sims)
    auto setSpeciesAmounts(VectorConstRef n) -> void;

    /// Set the species name to the new amount in xDC vector and adds this to bIC vector
    auto setSpeciesAmount(std::string name, double amount) -> void;

    /// Set the species with index ispecies to the new amount in xDC vector and adds this to bIC vector
    auto setSpeciesAmount(Index ispecies, double amount) -> void;

    /// Set the options of the ChemicalEngine instance
    auto setOptions(const ChemicalEngineOptions& options) -> void;

    /// Set smart start initial approximation (faster convergence, can be less accurate)
    auto setWarmStart() -> void;

    /// Set cold start automatic initial approximation (slower convergence, more accurate)
    auto setColdStart() -> void;

    // Caution: this may be ambiguous as in GEMS3K, species with the same name 
    //   may occur in more than one condensed phase!
    // An overload including phase name is needed!  
    /// Set the species @param name @param amount upper bound in dul vector
    /// if amount < 0 then resets to default 1e6
    auto setSpeciesUpperLimit(std::string name, double amount) -> void;

    /// Set the species @param name @param amount lower bound in dll vector
    /// if @param amount < 0 then resets to default 0
    auto setSpeciesLowerLimit(std::string name, double amount) -> void;

    /// Set species with index @param ispecies upper bound to @param amount in dul vector
    /// if @param amount < 0 then resets to default 1e6
    auto setSpeciesUpperLimit(Index ispecies, double amount) -> void;
   
    /// Set species with index @param ispecies upper bound to @param amount in dul vector
    /// if @param amount < 0 then resets to default 0 
    auto setSpeciesLowerLimit(Index ispecies, double amount) -> void;

    /// Set the species @param name StandardMolarGibbsEnergy @param value
    auto setStandardMolarGibbsEnergy(std::string name, double value) -> void;

    /// set all the upper limits for the species dul (in units of mol)
    auto setSpeciesUpperLimits(VectorConstRef n) -> void;
    /// set all the lower limits for the species dll (in units of mol)
    auto setSpeciesLowerLimits(VectorConstRef n) -> void;

  
    /// Sets pressure to @param P (in Pa) and temperature to @param T (in K)
    /// in ChemicalEngine instance without computing the equilibrium state
    /// Returns false if P.T was reset o.k., or true if P or T is out of 
    /// range for the lookup for thermodynamic data in the node (with old P,T retained).  
    auto setPT( double P, double T) const -> bool;  
    
    /// Sets the amounts of elements (vector b) in moles into ChemicalEngine
    /// instance without computing the equilibrium state
    auto setB( VectorConstRef b) -> void;

    /// recalculate the equilibrium state with preset P,T,B and other inputs
    /// If @param warmstart is true then GEMS3K uses contents of node as initial 
    ///   guess, otherwise (false) it gets the simplex LP initial guess.
    /// Returns the return code of GEMS3K GEM_run call
    auto reequilibrate(bool warmstart) -> int;

    /// Calculate the equilibrium state of the chemical system.
    /// @param T The temperature for the equilibrium calculation (in units of K)
    /// @param P The pressure for the equilibrium calculation (in units of Pa)
    /// @param b The amounts of the elements (in units of mol)
    /// Returns the return code of GEMS3K GEM_run call
    auto equilibrate(double T, double P, VectorConstRef b) -> int;

    /// Return the convergence result of the equilibrium calculation
    auto converged() const -> bool;

    /// Return the number of iterations of the equilibrium calculation
    auto numIterations() const -> Index;

    /// Return the wall time of the equilibrium calculation (in units of s)
    auto elapsedTime() const -> double;

    /// Return the temperature (in units of K)
    auto temperature() const -> double;

    /// Return the pressure (in units of Pa)
    auto pressure() const -> double;

    /// Return the amounts of the elements (in units of mol)
    auto elementAmounts() const -> VectorConstRef;

    /// Return the amounts of the elements in a phase.
    /// @param iphase The index of the phase.
    auto elementAmountsInPhase(Index iphase) const -> Vector;

    /// Return the amounts of the elements in a group of species.
    /// @param ispecies The vector of indices of the species.
    auto elementAmountsInSpecies(VectorXiConstRef ispecies) const -> Vector; // commented out until Eigen 3.4 is released on conda forge

    /// Return the amounts of the species (in units of mol)
    auto speciesAmounts() const -> VectorConstRef;

    /// Return the amount of the species with index @param ispecies (in units of mol)
    auto speciesAmount(Index ispecies) const -> double;

    /// Return the amount of the species with name @param name (in units of mol)
    auto speciesAmount(std::string name) const -> double;

    /// Return the upper limits for the species dul (in units of mol)
    auto speciesUpperLimits() const -> VectorConstRef;
    /// Return the lower limits for the species dll (in units of mol)
    auto speciesLowerLimits() const -> VectorConstRef;

    /// Return the molalities of the species. Aquatic systems only (assuming 
    /// the aqueous phase is the first one and H2O-solvent is the last species in it)
    auto speciesMolalities() const -> VectorConstRef;

    /// Return the mole fractions of the species.
    auto moleFractions() const -> VectorConstRef;

    /// Return the ln activity coefficients of the species (in mole fraction scale).
    auto lnActivityCoefficients() const -> VectorConstRef;

    /// Return the ln activities of the species.
    auto lnActivities() const -> VectorConstRef;

    /// Return the ln of concentrations of the species in their respective phases.
    auto lnConcentrations() const -> VectorConstRef;

    /// Return the chemical potentials of species (in units of J/mol).
    auto chemicalPotentials() const -> VectorConstRef;

    /// Return the standard  molar Gibbs energy of the species (in units of J/mol).
    auto standardMolarGibbsEnergy(Index ispecies) const -> double;

    /// Return the standard  molar enthalpy of the species (in units of J/mol).
    auto standardMolarEnthalpy(Index ispecies) const -> double;

    /// Return the standard  molar volume of the species (in units of m3/mol).
    auto standardMolarVolume(Index ispecies) const -> double;
    
    /// Return the standard  molar entropies of the species (in units of J/(mol*K)).
    auto standardMolarEntropy(Index ispecies) const -> double;

    /// TBD Return the standard  molar internal energies of the species (in units of J/mol).
    // auto standardMolarInternalEnergy(Index ispecies) const -> double;

    /// TBD Return the standard  molar Helmholtz energies of the species (in units of J/mol).
    // auto standardMolarHelmholtzEnergy(Index ispecies) const -> double;

    /// Return the standard  molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardMolarHeatCapacityConstP(Index ispecies) const -> double;

    /// TBD Return the standard  molar isochoric heat capacities of the species (in units of J/(mol*K)).
    // auto standardMolarHeatCapacityConstV(Index ispecies) const -> double;

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergy(Index iphase) const -> double;

    /// Return the molar enthalpy of the phase @param iphase (in units of J/mol).
    auto phaseMolarEnthalpy(Index iphase) const -> double;

    /// Return the molar volume of the phase @param iphase (in units of m3/mol).
    auto phaseMolarVolume(Index iphase) const -> double;

    /// Return the molar entropy of the phase @param iphase (in units of J/(mol*K)).
    auto phaseMolarEntropy(Index iphase) const -> double;

    /// TBD Return the molar internal energies of the phases (in units of J/mol).
    // auto phaseMolarInternalEnergy(Index iphase) const -> double;

    /// TBD Return the molar Helmholtz energies of the phases (in units of J/mol).
    // auto phaseMolarHelmholtzEnergy(Index iphase) const -> double;

    /// Return the molar isobaric heat capacity of phase @param iphase (in units of J/(mol*K)).
    auto phaseMolarHeatCapacityConstP(Index iphase) const -> double;

    /// TBD Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    // auto phaseMolarHeatCapacityConstV(Index ispecies) const -> double;

    /// Return the specific Gibbs energy of the phase @param iphase (in units of J/kg).
    auto phaseSpecificGibbsEnergy(Index iphase) const -> double;

    /// Return the specific enthalpy of the phase @param iphase (in units of J/kg).
    auto phaseSpecificEnthalpy(Index iphase) const -> double;

    /// Return the specific volume of the phase @param iphase (in units of m3/kg).
    auto phaseSpecificVolume(Index iphase) const -> double;

    /// Return the specific entropy of the phase @param iphase (in units of J/(kg*K)).
    auto phaseSpecificEntropy(Index iphase) const -> double;

    /// TBD Return the specific internal energies of the phases (in units of J/kg).
    // auto phaseSpecificInternalEnergy(Index iphase) const -> double;

    /// TBD Return the specific Helmholtz energies of the phases (in units of J/kg).
    // auto phaseSpecificHelmholtzEnergy(Index iphase) const -> double;

    /// Return the specific isobaric heat capacity of phase @param iphase (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacityConstP(Index iphase) const -> double;

    /// TBD Return the specific isochoric heat capacities of the phases (in units of J/(kg*K)).
    // auto phaseSpecificHeatCapacityConstV(Index iphase) const -> double;

    /// Return the densities of the phases (in units of kg/m3).
    auto phaseDensities() const -> VectorConstRef;

     /// Return the density of the phase @param iphase (in units of kg/m3).
    auto phaseDensity(Index iphase) const -> double;

    /// Return the masses of the phases (in units of kg).
    auto phaseMasses() const -> VectorConstRef;
 
    /// Return the mass of the phase @param iphase (in units of kg).
    auto phaseMass(Index iphase) const -> double;
    
    /// Return the molar amounts of the phases (in units of mol).
    auto phaseAmounts() const -> VectorConstRef;

    /// Return the molar amount of the phase @param iphase (in units of mol).
    auto phaseAmount(Index iphase) const -> double;

    /// Return the volumes of the phases (in units of m3).
    auto phaseVolumes() const -> VectorConstRef;
  
    /// Return the volume of the phase @param iphase (in units of m3).
    auto phaseVolume(Index iphase) const -> double;  

    /// Return the enthalpies of the phases (in units of J).
    auto phaseEnthalpies() const -> VectorConstRef;
    
    /// Returns the enthalpy of the phase @param iphase (in units of J).
    auto phaseEnthalpy(Index iphase) const -> double;

    /// Return the entropies of the phases (in units of J/K).
    auto phaseEntropies() const -> VectorConstRef;
    
    /// Returns the entropy of the phase @param iphase (in units of J/K).
    auto phaseEntropy(Index iphase) const -> double;

    /// Return the heat capacities Cp of the phases (in units of J/k).
    auto phaseHeatCapacitiesConstP() const -> VectorConstRef;
    
    /// Returns the heat capacity Cp of the phase @param iphase (in units of J/K).
    auto phaseHeatCapacityConstP(Index iphase) const -> double;

    /// Return the saturation (stability) indices  of the phases (in log10 units).
    auto phaseSatIndices() const -> VectorConstRef;
       
    /// Returns saturation (stability) index (in log10 units) of phase @param iphase.
    auto phaseSatIndex(Index iphase) const -> double;

    /// Return the mass of the system (in units of kg).
    auto systemMass() const -> double;

    /// Return the volume of the system (in units of m3).
    auto systemVolume() const -> double;

    /// Return the ionic strength of the aqueous phase (in units of molal).
    auto ionicStrength() const -> double;

    /// Return the pH of the aqueous phase.
    auto pH() const -> double;

    /// Return the pe of the aqueous phase.
    auto pe() const -> double;

    /// Return the Eh of the aqueous phase.
    auto Eh() const -> double;

    /// Return the total Gibbs energy of the system.
    auto systemGibbsEnergy() const -> double;

     /// Return the total enthalpy of the system.
    auto systemEnthalpy() const -> double;

    /// Return the total entropy of the system.
    auto systemEntropy() const -> double;

    /// Return the total heat capacity Cp of the system.
    auto systemHeatCapacityConstP() const -> double;

    /// The name of aqueous phase
    /// returns null string if no aq phase in system.
    auto aqueousPhaseName() const -> std::string;

    /// The name of gas phase
    /// returns null string if no gas phase in system.
    auto gasPhaseName() const -> std::string;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output the state of the ChemicalEngine object.
auto operator<<(std::ostream& out, const ChemicalEngine& engine) -> std::ostream&;

} // namespace xGEMS
