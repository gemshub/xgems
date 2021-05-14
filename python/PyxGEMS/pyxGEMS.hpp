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
        .def("initializeFromJsonStrings", &ChemicalEngine::initializeFromJsonStrings)
        .def("readDbrFromFile", &ChemicalEngine::readDbrFromFile)
        .def("readDbrFromJsonString", &ChemicalEngine::readDbrFromJsonString)
        .def("writeDbrToFile", &ChemicalEngine::writeDbrToFile)
        .def("writeDbrToJsonString", &ChemicalEngine::writeDbrToJsonString)
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
        .def("indexSpeciesAll", &ChemicalEngine::indexSpeciesAll)
        .def("indexPhase", &ChemicalEngine::indexPhase)
        .def("indexPhaseAll", &ChemicalEngine::indexPhaseAll)
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
        .def("setPT", &ChemicalEngine::setPT)
        .def("setB", &ChemicalEngine::setB)
        .def("reequilibrate", &ChemicalEngine::reequilibrate)
        .def("equilibrate", &ChemicalEngine::equilibrate)
        .def("converged", &ChemicalEngine::converged)
        .def("numIterations", &ChemicalEngine::numIterations)
        .def("elapsedTime", &ChemicalEngine::elapsedTime)
        .def("temperature", &ChemicalEngine::temperature)
        .def("pressure", &ChemicalEngine::pressure)
        .def("elementAmounts", &ChemicalEngine::elementAmounts, py::return_value_policy::reference_internal)
        .def("elementAmountsInPhase", &ChemicalEngine::elementAmountsInPhase)
        //.def("elementAmountsInSpecies", &ChemicalEngine::elementAmountsInSpecies)
        .def("speciesAmount", speciesAmount1)
        .def("speciesAmount", speciesAmount2)
        .def("setSpeciesAmounts", &ChemicalEngine::setSpeciesAmounts, py::return_value_policy::reference_internal)
        .def("speciesAmounts", &ChemicalEngine::speciesAmounts, py::return_value_policy::reference_internal)
        .def("setSpeciesUpperLimits", &ChemicalEngine::setSpeciesUpperLimits, py::return_value_policy::reference_internal)
        .def("setSpeciesLowerLimits", &ChemicalEngine::setSpeciesLowerLimits, py::return_value_policy::reference_internal)
        .def("speciesUpperLimits", &ChemicalEngine::speciesUpperLimits, py::return_value_policy::reference_internal)
        .def("speciesLowerLimits", &ChemicalEngine::speciesLowerLimits, py::return_value_policy::reference_internal)
        .def("moleFractions", &ChemicalEngine::moleFractions, py::return_value_policy::reference_internal)
        .def("speciesMolalities", &ChemicalEngine::speciesMolalities, py::return_value_policy::reference_internal)
        .def("lnActivityCoefficients", &ChemicalEngine::lnActivityCoefficients, py::return_value_policy::reference_internal)
        .def("lnActivities", &ChemicalEngine::lnActivities, py::return_value_policy::reference_internal)
        .def("lnConcentrations", &ChemicalEngine::lnConcentrations, py::return_value_policy::reference_internal)        
        .def("chemicalPotentials", &ChemicalEngine::chemicalPotentials, py::return_value_policy::reference_internal)
        .def("standardMolarGibbsEnergy", &ChemicalEngine::standardMolarGibbsEnergy)
        .def("standardMolarEnthalpy", &ChemicalEngine::standardMolarEnthalpy)
        .def("standardMolarVolume", &ChemicalEngine::standardMolarVolume)
        .def("standardMolarEntropy", &ChemicalEngine::standardMolarEntropy)
     //   .def("standardMolarInternalEnergy", &ChemicalEngine::standardMolarInternalEnergy)
     //   .def("standardMolarHelmholtzEnergy", &ChemicalEngine::standardMolarHelmholtzEnergy)
        .def("standardMolarHeatCapacityConstP", &ChemicalEngine::standardMolarHeatCapacityConstP)
     //   .def("standardMolarHeatCapacityConstV", &ChemicalEngine::standardMolarHeatCapacityConstV)
        .def("phaseMolarGibbsEnergy", &ChemicalEngine::phaseMolarGibbsEnergy)
        .def("phaseMolarEnthalpy", &ChemicalEngine::phaseMolarEnthalpy)
        .def("phaseMolarVolume", &ChemicalEngine::phaseMolarVolume)
        .def("phaseMolarEntropy", &ChemicalEngine::phaseMolarEntropy)
    //    .def("phaseMolarInternalEnergy", &ChemicalEngine::phaseMolarInternalEnergy)
    //    .def("phaseMolarHelmholtzEnergy", &ChemicalEngine::phaseMolarHelmholtzEnergy)
        .def("phaseMolarHeatCapacityConstP", &ChemicalEngine::phaseMolarHeatCapacityConstP)
    //    .def("phaseMolarHeatCapacityConstV", &ChemicalEngine::phaseMolarHeatCapacitiesConstV)
        .def("phaseSpecificGibbsEnergy", &ChemicalEngine::phaseSpecificGibbsEnergy)
        .def("phaseSpecificEnthalpy", &ChemicalEngine::phaseSpecificEnthalpy)
        .def("phaseSpecificVolume", &ChemicalEngine::phaseSpecificVolume)
        .def("phaseSpecificEntropy", &ChemicalEngine::phaseSpecificEntropy)
    //    .def("phaseSpecificInternalEnergy", &ChemicalEngine::phaseSpecificInternalEnergy)
    //    .def("phaseSpecificHelmholtzEnergy", &ChemicalEngine::phaseSpecificHelmholtzEnergy)
        .def("phaseSpecificHeatCapacityConstP", &ChemicalEngine::phaseSpecificHeatCapacityConstP)
    //    .def("phaseSpecificHeatCapacityConstV", &ChemicalEngine::phaseSpecificHeatCapacityConstVl)
        .def("phaseDensities", &ChemicalEngine::phaseDensities, py::return_value_policy::reference_internal)
        .def("phaseDensity", &ChemicalEngine::phaseDensity)
        .def("phaseMasses", &ChemicalEngine::phaseMasses, py::return_value_policy::reference_internal)
        .def("phaseMass", &ChemicalEngine::phaseMass)
        .def("phaseAmounts", &ChemicalEngine::phaseAmounts, py::return_value_policy::reference_internal)
        .def("phaseAmount", &ChemicalEngine::phaseAmount)
        .def("phaseVolumes", &ChemicalEngine::phaseVolumes, py::return_value_policy::reference_internal)
        .def("phaseVolume", &ChemicalEngine::phaseVolume)
        .def("phaseEnthalpies", &ChemicalEngine::phaseEnthalpies, py::return_value_policy::reference_internal)
        .def("phaseEnthalpy", &ChemicalEngine::phaseEnthalpy)
        .def("phaseEntropies", &ChemicalEngine::phaseEntropies, py::return_value_policy::reference_internal)
        .def("phaseEntropy", &ChemicalEngine::phaseEntropy)
        .def("phaseHeatCapacitiesConstP", &ChemicalEngine::phaseHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("phaseHeatCapacityConstP", &ChemicalEngine::phaseHeatCapacityConstP)
        .def("phaseSatIndices", &ChemicalEngine::phaseSatIndices, py::return_value_policy::reference_internal)
        .def("phaseSatIndex", &ChemicalEngine::phaseSatIndex)
        .def("systemMass", &ChemicalEngine::systemMass)
        .def("systemVolume", &ChemicalEngine::systemVolume)
        .def("systemGibbsEnergy", &ChemicalEngine::systemGibbsEnergy)
        .def("systemEnthalpy", &ChemicalEngine::systemEnthalpy)
        .def("systemEntropy", &ChemicalEngine::systemEntropy)
        .def("systemHeatCapacityConstP", &ChemicalEngine::systemHeatCapacityConstP)
        .def("ionicStrength", &ChemicalEngine::ionicStrength)
        .def("pH", &ChemicalEngine::pH)
        .def("pe", &ChemicalEngine::pe)
        .def("Eh", &ChemicalEngine::Eh)
        .def("__repr__", [](const ChemicalEngine& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
