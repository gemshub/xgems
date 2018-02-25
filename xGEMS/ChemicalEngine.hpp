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
#include <string>
#include <memory>

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
    /// Construct a default ChemicalEngine object
    ChemicalEngine();

    /// Construct a copy of a ChemicalEngine object
    ChemicalEngine(const ChemicalEngine& other);

    /// Construct a ChemicalEngine object from a specification file
    /// @param filename The name of the file containing the definition of the chemical system
    ChemicalEngine(std::string filename);

    /// Destroy this ChemicalEngine instance
    virtual ~ChemicalEngine();

    /// Assign another ChemicalEngine object to this
    auto operator=(ChemicalEngine other) -> ChemicalEngine&;

    /// Return the temperature (in units of K)
    auto temperature() const -> double;

    /// Return the pressure (in units of Pa)
    auto pressure() const -> double;

    /// Return the amounts of the elements (in units of mol)
    auto elementAmounts() const -> VectorConstRef;

    /// Return the amounts of the species (in units of mol)
    auto speciesAmounts() const -> VectorConstRef;

    /// Return the number of elements
    auto numElements() const -> unsigned;

    /// Return the number of species
    auto numSpecies() const -> unsigned;

    /// Return the number of phases
    auto numPhases() const -> unsigned;

    /// Return the number of species in a phase
    auto numSpeciesInPhase(Index iphase) const -> unsigned;

    /// Return the name of an element
    auto elementName(Index ielement) const -> std::string;

    /// Return the molar mass of an element (in units of kg/mol)
    auto elementMolarMass(Index ielement) const -> double;

    /// Return the stoichiometry of an element in a species
    auto elementStoichiometry(Index ispecies, Index ielement) const -> double;

    /// Return the name of a species
    auto speciesName(Index ispecies) const -> std::string;

    /// Return the name of a phase
    auto phaseName(Index iphase) const -> std::string;

    /// Set the options of the ChemicalEngine instance
    auto setOptions(const ChemicalEngineOptions& options) -> void;

    /// Calculate the equilibrium state of the system
    /// @param T The temperature for the equilibrium calculation (in units of K)
    /// @param P The pressure for the equilibrium calculation (in units of Pa)
    /// @param n The amounts of the elements (in units of mol)
    auto equilibrate(double T, double P, VectorConstRef b) -> void;

    /// Return the convergence result of the equilibrium calculation
    auto converged() const -> bool;

    /// Return the number of iterations of the equilibrium calculation
    auto numIterations() const -> unsigned;

    /// Return the wall time of the equilibrium calculation (in units of s)
    auto elapsedTime() const -> double;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace xGEMS
