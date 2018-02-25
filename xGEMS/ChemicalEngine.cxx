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
    py::class_<ChemicalEngine>(m, "ChemicalEngine")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("temperature", &ChemicalEngine::temperature)
        .def("pressure", &ChemicalEngine::pressure)
        .def("elementAmounts", &ChemicalEngine::elementAmounts, py::return_value_policy::reference_internal)
        .def("speciesAmounts", &ChemicalEngine::speciesAmounts, py::return_value_policy::reference_internal)
        .def("numElements", &ChemicalEngine::numElements)
        .def("numSpecies", &ChemicalEngine::numSpecies)
        .def("numPhases", &ChemicalEngine::numPhases)
        .def("numSpeciesInPhase", &ChemicalEngine::numSpeciesInPhase)
        .def("elementName", &ChemicalEngine::elementName)
        .def("elementMolarMass", &ChemicalEngine::elementMolarMass)
        .def("elementStoichiometry", &ChemicalEngine::elementStoichiometry)
        .def("speciesName", &ChemicalEngine::speciesName)
        .def("phaseName", &ChemicalEngine::phaseName)
        .def("setOptions", &ChemicalEngine::setOptions)
        .def("equilibrate", &ChemicalEngine::equilibrate)
        .def("converged", &ChemicalEngine::converged)
        .def("numIterations", &ChemicalEngine::numIterations)
        .def("elapsedTime", &ChemicalEngine::elapsedTime)
        ;
}
