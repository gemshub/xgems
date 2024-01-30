// ChemicalFun is a C++ and Python library for of C++ and Python API
// for Chemical Formula Parser and Reactions Generator.
//
// Copyright (C) 2018-2022 G.D.Miron, D.Kulik, S.Dmytriieva
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

#if _MSC_VER >= 1929
#include <corecrt.h>
#endif

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;
using namespace pybind11::literals;

// xGEMS includes
#include <xGEMS/GEMSEngine.h>
using namespace xGEMS;


void exportGEMSEngine(py::module& m)
{

    py::class_<GEMSEngine>(m, "GEMS2")
            .def(py::init<std::string, bool, bool>(), py::arg("input_file"), py::arg("reset_calc")=true, py::arg("cold_start")=true )
            .def_readwrite("T", &GEMSEngine::T)
            .def_readwrite("P", &GEMSEngine::P)
            .def_readwrite("b", &GEMSEngine::b_amounts)

            .def("element_names", &GEMSEngine::element_names)
            .def("species_names", &GEMSEngine::species_names)
            .def("phase_names", &GEMSEngine::phase_names)
            .def("aq_phase_symbol", &GEMSEngine::aq_phase_symbol)
            .def("gas_phase_symbol", &GEMSEngine::gas_phase_symbol)
            .def("element_molar_masses", &GEMSEngine::element_names)
            .def("species_in_phase", &GEMSEngine::species_in_phase)
            .def("species_charges", &GEMSEngine::species_charges)
            .def("species_molar_mass", &GEMSEngine::species_molar_mass)
            .def("species_molar_volumes", &GEMSEngine::species_molar_volumes)

            .def("equilibrate", &GEMSEngine::equilibrate)
            .def("cold_start", &GEMSEngine::cold_start)
            .def("warm_start", &GEMSEngine::warm_start)
            .def("nelements", &GEMSEngine::nelements)
            .def("nphases", &GEMSEngine::nphases)
            .def("nspecies", &GEMSEngine::nspecies)
            .def("clear", &GEMSEngine::clear, py::arg("cvalue")=1e-15)

            .def("set_species_G0", &GEMSEngine::set_species_G0)
            .def("bulk_composition", &GEMSEngine::bulk_composition)
            .def("pH", &GEMSEngine::pH)
            .def("pE", &GEMSEngine::pE)
            .def("ionic_strength", &GEMSEngine::ionic_strength)
            .def("system_volume", &GEMSEngine::system_volume)
            .def("system_mass", &GEMSEngine::system_mass)
            .def("phases_molar_volume", &GEMSEngine::phases_molar_volume)
            .def("phase_sat_indices", &GEMSEngine::phase_sat_indices)

            .def("aq_elements_molarity", &GEMSEngine::aq_elements_molarity)
            .def("aq_elements_molality", &GEMSEngine::aq_elements_molality)
            .def("aq_species_molarity", &GEMSEngine::aq_species_molarity)
            .def("aq_species_molality", &GEMSEngine::aq_species_molality)
            .def("aq_elements_moles", &GEMSEngine::aq_elements_moles)
            .def("set_bulk_composition", &GEMSEngine::set_bulk_composition)
            .def("reset_aq_composition", &GEMSEngine::reset_aq_composition)
            .def("solids_elements_moles", &GEMSEngine::solids_elements_moles)
            .def("phases_elements_moles", &GEMSEngine::phases_elements_moles)

            .def("phases_moles", &GEMSEngine::phases_moles)
            .def("species_moles", &GEMSEngine::species_moles)
            .def("species_ln_activities", &GEMSEngine::species_ln_activities)
            .def("species_ln_activity_coefficients", &GEMSEngine::species_ln_activity_coefficients)
            .def("phase_species_moles", &GEMSEngine::phase_species_moles)
            .def("solids_mass_frac", &GEMSEngine::solids_mass_frac)
            .def("solids_volume_frac", &GEMSEngine::solids_volume_frac)
            .def("aq_volume_frac", &GEMSEngine::aq_volume_frac)
            .def("phases_volume", &GEMSEngine::phases_volume)
            .def("phases_mass", &GEMSEngine::phases_mass)
            .def("phases_volume_frac", &GEMSEngine::phases_volume_frac)

            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_species_amt", &GEMSEngine::add_species_amt, py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("add_element_amt", &GEMSEngine::add_element_amt, py::arg("element_name"), py::arg("val"), py::arg("units")="moles")
            .def("add_multiple_elements_amt", &GEMSEngine::add_multiple_elements_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_amt_from_formula", &GEMSEngine::add_amt_from_formula, py::arg("formula"), py::arg("val"), py::arg("units")="moles")
            .def("get_b_from_formula", &GEMSEngine::get_b_from_formula, py::arg("formula"), py::arg("val")=1, py::arg("units")="moles")
            .def("set_multiple_species_lower_bound", &GEMSEngine::set_multiple_species_lower_bound, py::arg("input_dict"), py::arg("units")="moles")
            .def("set_multiple_species_upper_bound", &GEMSEngine::set_multiple_species_upper_bound, py::arg("input_dict"), py::arg("units")="moles")
            .def("set_species_lower_bound", &GEMSEngine::set_species_lower_bound, py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_upper_bound", &GEMSEngine::set_species_upper_bound, py::arg("species"), py::arg("val"), py::arg("units")="moles")

            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")

            .def("supress_phase", &GEMSEngine::supress_phase)
            .def("supress_multiple_phases", &GEMSEngine::supress_multiple_phases)
            .def("supress_species", &GEMSEngine::supress_species)
            .def("supress_multiple_species", &GEMSEngine::supress_multiple_species)
            .def("activate_phase", &GEMSEngine::activate_phase)
            .def("activate_multiple_phases", &GEMSEngine::activate_multiple_phases)
            .def("activate_species", &GEMSEngine::activate_species)
            .def("activate_multiple_species", &GEMSEngine::activate_multiple_species)
            ;

}
