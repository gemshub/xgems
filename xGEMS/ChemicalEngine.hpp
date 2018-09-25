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

    /// Initialize the ChemicalEngine object from a GEM-Selektor project file.
    /// @param filename The name of the file containing the definition of the chemical system
    auto initialize(std::string filename) -> void;

    /// Reads another DBR file (with input system composition, T,P etc.) \ . The DBR file must be compatible with
    /// the currently loaded IPM and DCH files (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFile(std::string filename) -> void;

    /// Return the number of elements.
    auto numElements() const -> Index;

    /// Return the number of species.
    auto numSpecies() const -> Index;

    /// Return the number of phases.
    auto numPhases() const -> Index;

    /// Return the number of species in a phase
    /// @param iphase The index of the phase.
    auto numSpeciesInPhase(Index iphase) const -> Index;

    /// Return the name of an element.
    /// @param ielement The index of the element.
    auto elementName(Index ielement) const -> std::string;

    /// Return the name of a species.
    /// @param ispecies The index of the species.
    auto speciesName(Index ispecies) const -> std::string;

    /// Return the name of a phase
    /// @param iphase The index of the phase.
    auto phaseName(Index iphase) const -> std::string;

    /// Return the index of an element.
    /// @param ielement The index of the element.
    auto indexElement(std::string element) const -> Index;

    /// Return the index of a species.
    /// @param ispecies The index of the species.
    auto indexSpecies(std::string species) const -> Index;

    /// Return the index of a phase
    /// @param iphase The index of the phase.
    auto indexPhase(std::string phase) const -> Index;

    /// Return the index of the phase with a species.
    /// @param ispecies The index of the species.
    auto indexPhaseWithSpecies(Index ispecies) const -> Index;

    /// Return the index of the first species in a phase.
    /// @param iphase The index of the phase.
    auto indexFirstSpeciesInPhase(Index iphase) const -> Index;

    /// Return the molar masses of the elements (in units of kg/mol)
    auto elementMolarMasses() const -> VectorConstRef;

    /// Return the molar masses of the species (in units of kg/mol)
    auto speciesMolarMasses() const -> VectorConstRef;

    /// Return the formula matrix of elements in species (in moles per mole species)
    /// (rows: elements, cols: species).
    auto formulaMatrix() const -> MatrixConstRef;

    /// Set the options of the ChemicalEngine instance
    auto setOptions(const ChemicalEngineOptions& options) -> void;

    /// Calculate the equilibrium state of the chemical system.
    /// @param T The temperature for the equilibrium calculation (in units of K)
    /// @param P The pressure for the equilibrium calculation (in units of Pa)
    /// @param b The amounts of the elements (in units of mol)
    auto equilibrate(double T, double P, VectorConstRef b) -> void;

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
    auto elementAmountsInSpecies(VectorXiConstRef ispecies) const -> Vector;

    /// Return the amounts of the species (in units of mol)
    auto speciesAmounts() const -> VectorConstRef;

    /// Return the molalities of the species.
    auto speciesMolalities() const -> VectorConstRef;

    /// Return the mole fractions of the species.
    auto moleFractions() const -> VectorConstRef;

    /// Return the ln activity coefficients of the species (in mole fraction scale).
    auto lnActivityCoefficients() const -> VectorConstRef;

    /// Return the ln activities of the species.
    auto lnActivities() const -> VectorConstRef;

    /// Return the ln of concentrations of the species in theri respective phases.
    auto lnConcentrations() const -> VectorConstRef;

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> VectorConstRef;

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> VectorConstRef;

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> VectorConstRef;

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> VectorConstRef;

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> VectorConstRef;

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> VectorConstRef;

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> VectorConstRef;

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> VectorConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> VectorConstRef;

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergies() const -> VectorConstRef;

    /// Return the molar enthalpies of the phases (in units of J/mol).
    auto phaseMolarEnthalpies() const -> VectorConstRef;

    /// Return the molar volumes of the phases (in units of m3/mol).
    auto phaseMolarVolumes() const -> VectorConstRef;

    /// Return the molar entropies of the phases (in units of J/(mol*K)).
    auto phaseMolarEntropies() const -> VectorConstRef;

    /// Return the molar internal energies of the phases (in units of J/mol).
    auto phaseMolarInternalEnergies() const -> VectorConstRef;

    /// Return the molar Helmholtz energies of the phases (in units of J/mol).
    auto phaseMolarHelmholtzEnergies() const -> VectorConstRef;

    /// Return the molar isobaric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstP() const -> VectorConstRef;

    /// Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstV() const -> VectorConstRef;

    /// Return the specific Gibbs energies of the phases (in units of J/kg).
    auto phaseSpecificGibbsEnergies() const -> VectorConstRef;

    /// Return the specific enthalpies of the phases (in units of J/kg).
    auto phaseSpecificEnthalpies() const -> VectorConstRef;

    /// Return the specific volumes of the phases (in units of m3/kg).
    auto phaseSpecificVolumes() const -> VectorConstRef;

    /// Return the specific entropies of the phases (in units of J/(kg*K)).
    auto phaseSpecificEntropies() const -> VectorConstRef;

    /// Return the specific internal energies of the phases (in units of J/kg).
    auto phaseSpecificInternalEnergies() const -> VectorConstRef;

    /// Return the specific Helmholtz energies of the phases (in units of J/kg).
    auto phaseSpecificHelmholtzEnergies() const -> VectorConstRef;

    /// Return the specific isobaric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstP() const -> VectorConstRef;

    /// Return the specific isochoric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstV() const -> VectorConstRef;

    /// Return the densities of the phases (in units of kg/m3).
    auto phaseDensities() const -> VectorConstRef;

    /// Return the masses of the phases (in units of kg).
    auto phaseMasses() const -> VectorConstRef;

    /// Return the molar amounts of the phases (in units of mol).
    auto phaseAmounts() const -> VectorConstRef;

    /// Return the volumes of the phases (in units of m3).
    auto phaseVolumes() const -> VectorConstRef;

    /// Return the saturation (stability) indices  of the phases (in log10 units).
    auto phaseSatIndices() const -> VectorConstRef;

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

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output the state of the ChemicalEngine object.
auto operator<<(std::ostream& out, const ChemicalEngine& engine) -> std::ostream&;

} // namespace xGEMS
