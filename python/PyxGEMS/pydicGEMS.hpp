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
namespace py=pybind11;
using namespace pybind11::literals;

// xGEMS includes
#include <xGEMS/GEMSEngine.hpp>
using namespace xGEMS;


void exportGEMSEngine(py::module& m)
{

    py::class_<GEMSEngine>  gems(m, "GEMS2");

    gems.def(py::init<std::string, bool, bool>(), py::arg("input_file"), py::arg("reset_calc")=true, py::arg("cold_start")=true )
            .def_readwrite("T", &GEMSEngine::T)
            .def_readwrite("P", &GEMSEngine::P)
            .def_readwrite("b", &GEMSEngine::b_amounts)

            .def_property_readonly("element_names", &GEMSEngine::element_names)
            .def_property_readonly("species_names", &GEMSEngine::species_names)
            .def_property_readonly("phase_names", &GEMSEngine::phase_names)
            .def_property_readonly("aq_phase_symbol", &GEMSEngine::aq_phase_symbol)
            .def_property_readonly("gas_phase_symbol", &GEMSEngine::gas_phase_symbol)
            .def_property_readonly("element_molar_masses", &GEMSEngine::element_molar_masses)
            .def_property_readonly("species_in_phase", &GEMSEngine::species_in_phase)
            .def_property_readonly("species_charges", &GEMSEngine::species_charges)
            .def_property_readonly("species_molar_mass", &GEMSEngine::species_molar_mass)
            .def_property_readonly("species_molar_volumes", &GEMSEngine::species_molar_volumes)

            .def("equilibrate", &GEMSEngine::equilibrate)
            .def("cold_start", &GEMSEngine::cold_start)
            .def("warm_start", &GEMSEngine::warm_start)
            .def_property_readonly("nelements", &GEMSEngine::nelements)
            .def_property_readonly("nphases", &GEMSEngine::nphases)
            .def_property_readonly("nspecies", &GEMSEngine::nspecies)
            .def("clear", &GEMSEngine::clear, py::arg("cvalue")=1e-15)

            .def("set_species_G0", &GEMSEngine::set_species_G0)
            .def_property_readonly("bulk_composition", &GEMSEngine::bulk_composition)
            //            .def_property_readonly("vector_b", &GEMSEngine::bulk_composition)
            .def_property_readonly("pH", &GEMSEngine::pH)
            .def_property_readonly("pE", &GEMSEngine::pE)
            .def_property_readonly("ionic_strength", &GEMSEngine::ionic_strength)
            .def_property_readonly("system_volume", &GEMSEngine::system_volume)
            .def_property_readonly("system_mass", &GEMSEngine::system_mass)
            .def_property_readonly("phases_molar_volume", &GEMSEngine::phases_molar_volume)
            .def_property_readonly("phase_sat_indices", &GEMSEngine::phase_sat_indices)

            .def_property_readonly("aq_elements_molarity", &GEMSEngine::aq_elements_molarity)
            .def_property_readonly("aq_elements_molality", &GEMSEngine::aq_elements_molality)
            .def_property_readonly("aq_species_molarity", &GEMSEngine::aq_species_molarity)
            .def_property_readonly("aq_species_molality", &GEMSEngine::aq_species_molality)
            .def_property_readonly("aq_elements_moles", &GEMSEngine::aq_elements_moles)
            .def("set_bulk_composition", &GEMSEngine::set_bulk_composition)
            .def("reset_aq_composition", &GEMSEngine::reset_aq_composition)
            .def_property_readonly("solids_elements_moles", &GEMSEngine::solids_elements_moles)
            .def_property_readonly("phases_elements_moles", &GEMSEngine::phases_elements_moles)

            .def("phase_species_moles", py::overload_cast<>(&GEMSEngine::phase_species_moles))
            .def("phase_species_moles", py::overload_cast<std::string>(&GEMSEngine::phase_species_moles))
            .def_property_readonly("phases_moles", &GEMSEngine::phases_moles)
            .def_property_readonly("species_moles", &GEMSEngine::species_moles)
            .def_property_readonly("species_ln_activities", &GEMSEngine::species_ln_activities)
            .def_property_readonly("species_ln_activity_coefficients", &GEMSEngine::species_ln_activity_coefficients)
            .def_property_readonly("species_upper_bounds", &GEMSEngine::species_upper_bounds)
            .def_property_readonly("species_lower_bounds", &GEMSEngine::species_lower_bounds)
            .def_property_readonly("phase_species_ln_activities", &GEMSEngine::phase_species_ln_activities)
            .def_property_readonly("phase_species_ln_activity_coefficients", &GEMSEngine::phase_species_ln_activity_coefficients)
            .def_property_readonly("phase_species_upper_bounds", &GEMSEngine::phase_species_upper_bounds)
            .def_property_readonly("phase_species_lower_bounds", &GEMSEngine::phase_species_lower_bounds)
            .def_property_readonly("solids_mass_frac", &GEMSEngine::solids_mass_frac)
            .def_property_readonly("solids_volume_frac", &GEMSEngine::solids_volume_frac)
            .def_property_readonly("aq_volume_frac", &GEMSEngine::aq_volume_frac)
            .def_property_readonly("phases_volume", &GEMSEngine::phases_volume)
            .def_property_readonly("phases_mass", &GEMSEngine::phases_mass)
            .def_property_readonly("phases_volume_frac", &GEMSEngine::phases_volume_frac)

            .def("add_multiple_species_amt", &GEMSEngine::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_species_amt", &GEMSEngine::add_species_amt, py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("add_element_amt", &GEMSEngine::add_element_amt, py::arg("element_name"), py::arg("val"), py::arg("units")="moles")
            .def("add_multiple_elements_amt", &GEMSEngine::add_multiple_elements_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_amt_from_formula", &GEMSEngine::add_amt_from_formula, py::arg("formula"), py::arg("val"), py::arg("units")="moles")
            .def("get_b_from_formula", &GEMSEngine::get_b_from_formula, py::arg("formula"), py::arg("val")=1, py::arg("units")="moles")

            .def("set_multiple_species_lower_bound", &GEMSEngine::set_multiple_species_lower_bound, py::arg("input_dict"), py::arg("units")="moles")
            .def("set_multiple_species_upper_bound", &GEMSEngine::set_multiple_species_upper_bound, py::arg("input_dict"), py::arg("units")="moles")
            .def("set_species_lower_bound", py::overload_cast<const std::string&, double, const std::string&>
                 (&GEMSEngine::set_species_lower_bound), py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_lower_bound", py::overload_cast<Index, double, const std::string&>
                 (&GEMSEngine::set_species_lower_bound), py::arg("ispecies"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_upper_bound", py::overload_cast<const std::string&, double, const std::string&>
                 (&GEMSEngine::set_species_upper_bound), py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_upper_bound", py::overload_cast<Index, double, const std::string&>
                 (&GEMSEngine::set_species_upper_bound), py::arg("ispecies"), py::arg("val"), py::arg("units")="moles")

            .def("supress_phase", &GEMSEngine::supress_phase)
            .def("supress_multiple_phases", &GEMSEngine::supress_multiple_phases)
            .def("supress_species", &GEMSEngine::supress_species)
            .def("supress_multiple_species", &GEMSEngine::supress_multiple_species)
            .def("activate_phase", &GEMSEngine::activate_phase)
            .def("activate_multiple_phases", &GEMSEngine::activate_multiple_phases)
            .def("activate_species", &GEMSEngine::activate_species)
            .def("activate_multiple_species", &GEMSEngine::activate_multiple_species)
            ;

    gems.attr("vector_b")  =   gems.attr("bulk_composition");
    gems.attr("aq_el_M")  = gems.attr("aq_elements_molarity");
    gems.attr("aq_el_my") = gems.attr("aq_elements_molality");
    gems.attr("aq_species_composition") = gems.attr("aq_elements_molality");  // this alias is misleading!
    gems.attr("aq_sp_M") = gems.attr("aq_species_molarity");
    gems.attr("aq_sp_my") = gems.attr("aq_species_molality");
    gems.attr("aq_composition") = gems.attr("aq_species_molality");   // this alias is misleading!
    gems.attr("aq_elements_amounts") = gems.attr("aq_elements_moles");
    gems.attr("set_vector_b") = gems.attr("set_bulk_composition");
    gems.attr("clear_b_aq_part") = gems.attr("reset_aq_composition");
    gems.attr("solid_elements_amounts") = gems.attr("solids_elements_moles");
    gems.attr("phase_elements_amounts") = gems.attr("phases_elements_moles");
    gems.attr("phase_amounts") = gems.attr("phases_moles");
    gems.attr("species_amounts") = gems.attr("species_moles");
    gems.attr("phase_species_amounts") = gems.attr("phase_species_moles");
    gems.attr("solid_mass_frac") = gems.attr("solids_mass_frac");
    gems.attr("solid_volume_frac") = gems.attr("solids_volume_frac");
    gems.attr("phase_volumes") = gems.attr("phases_volume");
    gems.attr("phase_masses") = gems.attr("phases_mass");
    gems.attr("phase_volume_frac") = gems.attr("phases_volume_frac");
    gems.attr("vector_b_from_formula") = gems.attr("get_b_from_formula");
    gems.attr("multiple_species_lower_bound") = gems.attr("set_multiple_species_lower_bound");
    gems.attr("multiple_species_upper_bound") = gems.attr("set_multiple_species_upper_bound");
    gems.attr("IS") = gems.attr("ionic_strength");
    gems.attr("phase_molar_volume") = gems.attr("phases_molar_volume");
    gems.attr("phases_sat_index") = gems.attr("phase_sat_indices");
    gems.attr("phases_SI") = gems.attr("phase_sat_indices");

}
