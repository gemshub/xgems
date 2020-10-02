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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// xGEMS includes
#include <xGEMS/ChemicalEngine.hpp>
using namespace xGEMS;

void exportChemicalEngine(py::module& m)
{

    auto speciesAmount1 = static_cast<double (ChemicalEngine::*)(Index) const> (&ChemicalEngine::speciesAmount);
    auto speciesAmount2 = static_cast<double (ChemicalEngine::*)(std::string) const> (&ChemicalEngine::speciesAmount);

    auto setSpeciesAmount1 = static_cast<void (ChemicalEngine::*)(std::string, double)> (&ChemicalEngine::setSpeciesAmount);
    auto setSpeciesAmount2 = static_cast<void (ChemicalEngine::*)(Index, double)> (&ChemicalEngine::setSpeciesAmount);

    auto setSpeciesUpperLimit1 = static_cast<void (ChemicalEngine::*)(std::string, double)> (&ChemicalEngine::setSpeciesUpperLimit);
    auto setSpeciesUpperLimit2 = static_cast<void (ChemicalEngine::*)(Index, double)> (&ChemicalEngine::setSpeciesUpperLimit);

    auto setSpeciesLowerLimit1 = static_cast<void (ChemicalEngine::*)(std::string, double)> (&ChemicalEngine::setSpeciesLowerLimit);
    auto setSpeciesLowerLimit2 = static_cast<void (ChemicalEngine::*)(Index, double)> (&ChemicalEngine::setSpeciesLowerLimit);

    py::class_<ChemicalEngine>(m, "ChemicalEngine")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("initialize", &ChemicalEngine::initialize)
        .def("readDbrFile", &ChemicalEngine::readDbrFile)
        .def("numElements", &ChemicalEngine::numElements)
        .def("numSpecies", &ChemicalEngine::numSpecies)
        .def("numPhases", &ChemicalEngine::numPhases)
        .def("numSpeciesInPhase", &ChemicalEngine::numSpeciesInPhase)
        .def("elementName", &ChemicalEngine::elementName)
        .def("speciesName", &ChemicalEngine::speciesName)
        .def("speciesCharge", &ChemicalEngine::speciesCharge)
        .def("phaseName", &ChemicalEngine::phaseName)
        .def("indexElement", &ChemicalEngine::indexElement)
        .def("indexSpecies", &ChemicalEngine::indexSpecies)
        .def("indexPhase", &ChemicalEngine::indexPhase)
        .def("indexPhaseWithSpecies", &ChemicalEngine::indexPhaseWithSpecies)
        .def("indexFirstSpeciesInPhase", &ChemicalEngine::indexFirstSpeciesInPhase)
        .def("elementMolarMasses", &ChemicalEngine::elementMolarMasses, py::return_value_policy::reference_internal)
        .def("speciesMolarMasses", &ChemicalEngine::speciesMolarMasses, py::return_value_policy::reference_internal)
        .def("formulaMatrix", &ChemicalEngine::formulaMatrix, py::return_value_policy::reference_internal)
        .def("setOptions", &ChemicalEngine::setOptions)
        .def("setWarmStart", &ChemicalEngine::setWarmStart)
        .def("setColdStart", &ChemicalEngine::setColdStart)
        .def("setSpeciesUpperLimit", setSpeciesUpperLimit1)
        .def("setSpeciesUpperLimit", setSpeciesUpperLimit2)
        .def("setSpeciesLowerLimit", setSpeciesLowerLimit1)
        .def("setSpeciesLowerLimit", setSpeciesLowerLimit2)
        .def("setSpeciesAmount", setSpeciesAmount1)
        .def("setSpeciesAmount", setSpeciesAmount2)
        .def("equilibrate", &ChemicalEngine::equilibrate)
        .def("converged", &ChemicalEngine::converged)
        .def("numIterations", &ChemicalEngine::numIterations)
        .def("elapsedTime", &ChemicalEngine::elapsedTime)
        .def("temperature", &ChemicalEngine::temperature)
        .def("pressure", &ChemicalEngine::pressure)
        .def("elementAmounts", &ChemicalEngine::elementAmounts, py::return_value_policy::reference_internal)
        .def("elementAmountsInPhase", &ChemicalEngine::elementAmountsInPhase)
        .def("elementAmountsInSpecies", &ChemicalEngine::elementAmountsInSpecies)
        .def("speciesAmount", speciesAmount1)
        .def("speciesAmount", speciesAmount2)
        .def("setSpeciesAmounts", &ChemicalEngine::setSpeciesAmounts, py::return_value_policy::reference_internal)
        .def("speciesAmounts", &ChemicalEngine::speciesAmounts, py::return_value_policy::reference_internal)
        .def("moleFractions", &ChemicalEngine::moleFractions, py::return_value_policy::reference_internal)
        .def("speciesMolalities", &ChemicalEngine::speciesMolalities, py::return_value_policy::reference_internal)
        .def("lnActivityCoefficients", &ChemicalEngine::lnActivityCoefficients, py::return_value_policy::reference_internal)
        .def("lnActivities", &ChemicalEngine::lnActivities, py::return_value_policy::reference_internal)
        .def("lnConcentrations", &ChemicalEngine::lnConcentrations, py::return_value_policy::reference_internal)        
        .def("chemicalPotentials", &ChemicalEngine::chemicalPotentials, py::return_value_policy::reference_internal)
        .def("standardMolarGibbsEnergies", &ChemicalEngine::standardMolarGibbsEnergies, py::return_value_policy::reference_internal)
        .def("standardMolarEnthalpies", &ChemicalEngine::standardMolarEnthalpies, py::return_value_policy::reference_internal)
        .def("standardMolarVolumes", &ChemicalEngine::standardMolarVolumes, py::return_value_policy::reference_internal)
        .def("standardMolarEntropies", &ChemicalEngine::standardMolarEntropies, py::return_value_policy::reference_internal)
        .def("standardMolarInternalEnergies", &ChemicalEngine::standardMolarInternalEnergies, py::return_value_policy::reference_internal)
        .def("standardMolarHelmholtzEnergies", &ChemicalEngine::standardMolarHelmholtzEnergies, py::return_value_policy::reference_internal)
        .def("standardMolarHeatCapacitiesConstP", &ChemicalEngine::standardMolarHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("standardMolarHeatCapacitiesConstV", &ChemicalEngine::standardMolarHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        .def("phaseMolarGibbsEnergies", &ChemicalEngine::phaseMolarGibbsEnergies, py::return_value_policy::reference_internal)
        .def("phaseMolarEnthalpies", &ChemicalEngine::phaseMolarEnthalpies, py::return_value_policy::reference_internal)
        .def("phaseMolarVolumes", &ChemicalEngine::phaseMolarVolumes, py::return_value_policy::reference_internal)
        .def("phaseMolarEntropies", &ChemicalEngine::phaseMolarEntropies, py::return_value_policy::reference_internal)
        .def("phaseMolarInternalEnergies", &ChemicalEngine::phaseMolarInternalEnergies, py::return_value_policy::reference_internal)
        .def("phaseMolarHelmholtzEnergies", &ChemicalEngine::phaseMolarHelmholtzEnergies, py::return_value_policy::reference_internal)
        .def("phaseMolarHeatCapacitiesConstP", &ChemicalEngine::phaseMolarHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("phaseMolarHeatCapacitiesConstV", &ChemicalEngine::phaseMolarHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        .def("phaseSpecificGibbsEnergies", &ChemicalEngine::phaseSpecificGibbsEnergies, py::return_value_policy::reference_internal)
        .def("phaseSpecificEnthalpies", &ChemicalEngine::phaseSpecificEnthalpies, py::return_value_policy::reference_internal)
        .def("phaseSpecificVolumes", &ChemicalEngine::phaseSpecificVolumes, py::return_value_policy::reference_internal)
        .def("phaseSpecificEntropies", &ChemicalEngine::phaseSpecificEntropies, py::return_value_policy::reference_internal)
        .def("phaseSpecificInternalEnergies", &ChemicalEngine::phaseSpecificInternalEnergies, py::return_value_policy::reference_internal)
        .def("phaseSpecificHelmholtzEnergies", &ChemicalEngine::phaseSpecificHelmholtzEnergies, py::return_value_policy::reference_internal)
        .def("phaseSpecificHeatCapacitiesConstP", &ChemicalEngine::phaseSpecificHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("phaseSpecificHeatCapacitiesConstV", &ChemicalEngine::phaseSpecificHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        .def("phaseDensities", &ChemicalEngine::phaseDensities, py::return_value_policy::reference_internal)
        .def("phaseMasses", &ChemicalEngine::phaseMasses, py::return_value_policy::reference_internal)
        .def("phaseAmounts", &ChemicalEngine::phaseAmounts, py::return_value_policy::reference_internal)
        .def("phaseVolumes", &ChemicalEngine::phaseVolumes, py::return_value_policy::reference_internal)
        .def("phaseVolume", &ChemicalEngine::phaseVolume, py::return_value_policy::reference_internal)
        .def("phaseMass", &ChemicalEngine::phaseMass, py::return_value_policy::reference_internal)
        .def("phaseEnthalpies", &ChemicalEngine::phaseEnthalpies, py::return_value_policy::reference_internal)
        .def("phaseEnthalpy", &ChemicalEngine::phaseEnthalpy, py::return_value_policy::reference_internal)
        .def("phaseSatIndices", &ChemicalEngine::phaseSatIndices, py::return_value_policy::reference_internal)
        .def("systemMass", &ChemicalEngine::systemMass)
        .def("systemVolume", &ChemicalEngine::systemVolume)
        .def("ionicStrength", &ChemicalEngine::ionicStrength)
        .def("pH", &ChemicalEngine::pH)
        .def("pe", &ChemicalEngine::pe)
        .def("Eh", &ChemicalEngine::Eh)
        .def("systemGibbsEnergy", &ChemicalEngine::systemGibbsEnergy)
        .def("systemEnthalpy", &ChemicalEngine::systemEnthalpy)
        .def("systemEntropy", &ChemicalEngine::systemEntropy)
        .def("systemHeatCapacityConstP", &ChemicalEngine::systemHeatCapacityConstP)
        .def("__repr__", [](const ChemicalEngine& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
