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
#include <xGEMS/ChemicalEngineMaps.hpp>
using namespace xGEMS;


void exportChemicalEngineMaps(py::module& m)
{

    py::class_<ChemicalEngineMaps>  gems(m, "ChemicalEngineDicts",
                                        R"doc(
Class for equilibrium computations and thermodynamic analysis using dictionaries.
Gems interface in calculator format for easy  using dictionaries.)doc");

    gems.def(py::init<std::string, bool, bool>(), py::arg("input_file"), py::arg("reset_calc")=false, py::arg("cold_start")=true,
             R"doc(
  Constructs a ChemicalEngineDicts instance by loading a GEM-Selektor project file.

  :param str input_file: The file path for the chemical system definition (e.g., "my-system-dat.lst").
  :param bool reset_calc: If true, clear the amounts of all elements, default false.
  :param bool cold_start: If true, configures the engine to use a cold start, default true.

  **Example:**

  .. code-block:: python

      engine = ChemicalEngineDicts("my-system-dat.lst")
  )doc" )
            .def_readwrite("T", &ChemicalEngineMaps::T,
             R"doc(
   Sets and gets the temperature without computing equilibrium.
   Property temperature in Kelvin (K).

  **Example:**

  .. code-block:: python

      engine.T = 298.15
  )doc")
            .def_readwrite("P", &ChemicalEngineMaps::P,
                            R"doc(
    Sets and gets the pressure without computing equilibrium.
    Property pressure in Pascals (Pa)

  **Example:**

  .. code-block:: python

      engine.P = 100000.0
  )doc")
            .def_property_readonly("element_names", &ChemicalEngineMaps::element_names,
                               R"doc(
    Read-only property: the names of all elements in the system.

  **Example:**

  .. code-block:: python

      print("element_names", engine.element_names)
  )doc")
            .def_property_readonly("species_names", &ChemicalEngineMaps::species_names,
                               R"doc(
    Read-only property: the names of all species in the system.

  **Example:**

  .. code-block:: python

      print("species_names", engine.species_names)
  )doc")
            .def_property_readonly("phase_names", &ChemicalEngineMaps::phase_names,
                               R"doc(
    Read-only property: the names of all phases in the system.

  **Example:**

  .. code-block:: python

      print("species_names", engine.species_names)
  )doc")
            .def_property_readonly("aq_phase_symbol", &ChemicalEngineMaps::aq_phase_symbol,
                               R"doc(
    Read-only property: the aqueous phase name. If empty, the aqueous phase is not in system.

  **Example:**

  .. code-block:: python

      print("aq_phase_symbol", engine.aq_phase_symbol)
  )doc")
            .def_property_readonly("gas_phase_symbol", &ChemicalEngineMaps::gas_phase_symbol,
                               R"doc(
    Read-only property: the gaseous phase name. If empty, the gaseous phase is not in system.

  **Example:**

  .. code-block:: python

      print("gas_phase_symbol", engine.gas_phase_symbol)
  )doc")
            .def_property_readonly("element_molar_masses", &ChemicalEngineMaps::element_molar_masses,
                               R"doc(
    Read-only property: the dictionary of molar masses (kg/mol) for each element.

  **Example:**

  .. code-block:: python

      print("element_molar_masses", engine.element_molar_masses)
  )doc")
            .def_property_readonly("species_in_phase", &ChemicalEngineMaps::species_in_phase,
                               R"doc(
    Read-only property: the dictionary of the names of all species for each phase in the system.

  **Example:**

  .. code-block:: python

      out = engine.species_in_phase
      for name in engine.phase_names:  #print fixed order to compare
          print(name, out[name])
  )doc")
            .def_property_readonly("species_charges", &ChemicalEngineMaps::species_charges,
                               R"doc(
    Read-only property: the dictionary of the electrical charge of a species.

  **Example:**

  .. code-block:: python

      out = engine.species_charges
      for name in engine.species_names:  #print fixed order to compare
          print(name, out[name])
  )doc")
            .def_property_readonly("species_molar_mass", &ChemicalEngineMaps::species_molar_mass,
                               R"doc(
    Read-only property: the dictionary of the species molar masses (kg/mol).

  **Example:**

  .. code-block:: python

      out = engine.species_molar_mass
      for name in engine.species_names:  #print fixed order to compare
          print(name, out[name])
  )doc")
            .def_property_readonly("species_molar_volumes", &ChemicalEngineMaps::species_molar_volumes,
                               R"doc(
    Read-only property: the dictionary of the species standard molar volumes in m³/mol.

  **Example:**

  .. code-block:: python

      out = engine.species_molar_volumes
      for name in engine.species_names:  #print fixed order to compare
          print(name, out[name])
  )doc")

        .def("equilibrate", &ChemicalEngineMaps::equilibrate,
             R"doc(
    Computes the equilibrium stateof the current system.
    Uses current temperature (K), pressure (Pa), and element amounts (in mol) to compute equilibrium.

  **Example:**

  .. code-block:: python

      engine.T = 298.15
      engine.P = 100000.0
      bulk_composition = {'C': 1e-08, 'Ca': 1e-08, 'Cl': 0.002, 'H': 111.016746657646,
                          'Mg': 0.001, 'O': 55.5083933588231, 'Sn': 130.841288437146, 'Zz': 0.0}
      engine.set_bulk_composition(bulk_composition)
      engine.equilibrate();

  **Return**
     The function returns the string indicating the status:

    - No GEM re-calculation needed
    - Need GEM calculation with LPP (automatic) initial approximation (AIA)
    - OK after GEM calculation with LPP AIA
    - Bad (not fully trustful) result after GEM calculation with LPP AIA
    - Failure (no result) in GEM calculation with LPP AIA
    - Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
    - OK after GEM calculation with SIA
    - Bad (not fully trustful) result after GEM calculation with SIA
    - Failure (no result) in GEM calculation with SIA
    - Terminal error in GEMS3K (e.g., memory corruption). Restart required.
  )doc")
            .def("cold_start", &ChemicalEngineMaps::cold_start,
             R"doc(
Enables cold start for the ChemicalEngineDicts.

Cold start resets the initial guess to default values for robust convergence.

**Example:**

.. code-block:: python

  engine.cold_start()
)doc")
            .def("warm_start", &ChemicalEngineMaps::warm_start,
             R"doc(
Enables warm start for the ChemicalEngineDicts.

Warm start uses the previous solution as an initial guess for faster convergence.

**Example:**

.. code-block:: python

  engine.warm_start()
)doc")
            .def_property_readonly("nelements", &ChemicalEngineMaps::nelements,
                               R"doc(
    Read-only property: the number of elements in the system.

  **Example:**

  .. code-block:: python

      print("nelements", engine.nelements)
  )doc")
            .def_property_readonly("nphases", &ChemicalEngineMaps::nphases,
                               R"doc(
    Read-only property: the number of phases in the system.

  **Example:**

  .. code-block:: python

      print("nphases", engine.nphases)
  )doc")
            .def_property_readonly("nspecies", &ChemicalEngineMaps::nspecies,
                               R"doc(
    Read-only property: the number of species in the system.

  **Example:**

  .. code-block:: python

      print("nspecies", engine.nspecies)
  )doc")
            .def("clear", &ChemicalEngineMaps::clear, py::arg("min_amount")=1e-15,
             R"doc(
Clear the amounts of elements (set the default amount for all components).

:param float min_amount: min amount of element in mole, default 1e-15.

**Example:**

.. code-block:: python

  engine.clear(1e-12)
)doc")

            .def("set_species_G0", &ChemicalEngineMaps::set_species_G0,
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).

:param str name: Species name.
:param float value: Standard molar Gibbs energy value (J/mol).

**Example:**

.. code-block:: python

  engine.set_species_G0("H2O", -237.13)
)doc")
            .def_property_readonly("bulk_composition", &ChemicalEngineMaps::bulk_composition,
                               R"doc(
    Read-only property: the dictionary of the amounts of the elements in mol.

  **Example:**

  .. code-block:: python

      print("bulk_composition", engine.bulk_composition)
  )doc")
            .def_property_readonly("pH", &ChemicalEngineMaps::pH,
                               R"doc(
    Read-only property: the pH of the aqueous phase (in the activity scale (-log10 molal)).

  **Example:**

  .. code-block:: python

      print("pH", engine.pH)
  )doc")
            .def_property_readonly("pE", &ChemicalEngineMaps::pE,
                               R"doc(
    Read-only property: the pe of the aqueous phase (in the activity scale (-log10 molal)).

  **Example:**

  .. code-block:: python

      print("pE", engine.pE)
  )doc")
            .def_property_readonly("ionic_strength", &ChemicalEngineMaps::ionic_strength,
                               R"doc(
    Read-only property: the ionic strength of the aqueous phase in molal.

  **Example:**

  .. code-block:: python

      print("ionic_strength", engine.ionic_strength)
  )doc")
            .def_property_readonly("system_volume", &ChemicalEngineMaps::system_volume,
                               R"doc(
    Read-only property: the total volume of the system in m³.

  **Example:**

  .. code-block:: python

      print("system_volume", engine.system_volume)
  )doc")
            .def_property_readonly("system_mass", &ChemicalEngineMaps::system_mass,
                               R"doc(
    Read-only property: the total mass of the system in kg.

  **Example:**

  .. code-block:: python

      print("system_mass", engine.system_mass)
  )doc")
            .def_property_readonly("phases_molar_volume", &ChemicalEngineMaps::phases_molar_volume,
                               R"doc(
    Read-only property: the dictionary of the phases molar volumes in m³/mol.

  **Example:**

  .. code-block:: python

      print("phases_molar_volume", engine.phases_molar_volume)
  )doc")
            .def_property_readonly("phase_sat_indices", &ChemicalEngineMaps::phase_sat_indices,
                               R"doc(
    Read-only property: the dictionary of the saturation indices of all phases (log₁₀ units).

  **Example:**

  .. code-block:: python

      print("phase_sat_indices", engine.phase_sat_indices)
  )doc")

            .def_property_readonly("aq_elements_molarity", &ChemicalEngineMaps::aq_elements_molarity,
                               R"doc(
    Read-only property: the dictionary for aq elements the aq solution composition in mol/L aq solution.

  **Example:**

  .. code-block:: python

      print("aq_elements_molarity", engine.aq_elements_molarity)
  )doc")
            .def_property_readonly("aq_elements_molality", &ChemicalEngineMaps::aq_elements_molality,
                               R"doc(
    Read-only property: the dictionary for aq elements the aq solution elemental composition in mol/kgH2O.

  **Example:**

  .. code-block:: python

      print("aq_elements_molality", engine.aq_elements_molality)
  )doc")
            .def_property_readonly("aq_species_molarity", &ChemicalEngineMaps::aq_species_molarity,
                               R"doc(
    Read-only property: the dictionary for aq species the aq solution composition in mol/L of aqueous species.

  **Example:**

  .. code-block:: python

      print("aq_species_molarity", engine.aq_species_molarity)
  )doc")
            .def_property_readonly("aq_species_molality", &ChemicalEngineMaps::aq_species_molality,
                               R"doc(
    Read-only property: the dictionary for the aq solution composition in mol/kg H2O of aqueous species (speciation).

  **Example:**

  .. code-block:: python

      print("aq_species_molality", engine.aq_species_molality)
  )doc")
            .def_property_readonly("aq_elements_moles", &ChemicalEngineMaps::aq_elements_moles,
                               R"doc(
    Read-only property: the dictionary of the amounts of each element in the aqueous phase (in mol).

  **Example:**

  .. code-block:: python

      print("aq_elements_moles", engine.aq_elements_moles)
  )doc")
            .def("set_bulk_composition", &ChemicalEngineMaps::set_bulk_composition, py::arg("b_input"), py::arg("min_amount")=1e-15,
             R"doc(
Sets the amounts of elements (vector b).

:param dict b_input: dictionary of elements amounts in mol.
:param float min_amount: min amount of element in mol, default 1e-15

**Example:**

.. code-block:: python

  bulk_composition = {'C': 1e-08, 'Ca': 1e-08, 'Cl': 0.002, 'H': 111.016746657646,
                    'Mg': 0.001, 'O': 55.5083933588231, 'Sn': 130.841288437146, 'Zz': 0.0}
  engine.set_bulk_composition(bulk_composition)
)doc")
            .def("reset_aq_composition", &ChemicalEngineMaps::reset_aq_composition, py::arg("min_amount")=1e-15,
             R"doc(
Removes bulk elemental aqueous solution composition from vector b.
Be careful as this will also remove water i.e H+ and OH-.

:param float min_amount: min amount of element in mol, default 1e-15

**Example:**

.. code-block:: python

  engine.reset_aq_composition()
)doc")
            .def("solids_elements_moles", &ChemicalEngineMaps::solids_elements_moles, py::arg("min_amount_phase")=1e-12, py::arg("min_amount_element")=1e-14,
                               R"doc(
    Returns the dictionary of the mole amounts of elements in all solids together.

  :param float min_amount_phase: min amount of phase in mol, default 1e-12.
  :param float min_amount_element: min amount of element in mol, default 1e-14.

  **Example:**

  .. code-block:: python

      print("solids_elements_moles", engine.solids_elements_moles())
  )doc")
            .def_property_readonly("phases_elements_moles", &ChemicalEngineMaps::phases_elements_moles,
                               R"doc(
    Read-only property: the dictionary of dictionaries containing mole amounts of elements for each phase (in mol).

  **Example:**

  .. code-block:: python

      print("phases_elements_moles", engine.phases_elements_moles)
  )doc")

            .def("phase_species_moles", py::overload_cast<>(&ChemicalEngineMaps::phase_species_moles))
            .def("phase_species_moles", py::overload_cast<std::string>(&ChemicalEngineMaps::phase_species_moles))
            .def_property_readonly("phases_moles", &ChemicalEngineMaps::phases_moles)
            .def_property_readonly("species_moles", &ChemicalEngineMaps::species_moles)
            .def_property_readonly("species_ln_activities", &ChemicalEngineMaps::species_ln_activities)
            .def_property_readonly("species_ln_activity_coefficients", &ChemicalEngineMaps::species_ln_activity_coefficients)
            .def_property_readonly("species_upper_bounds", &ChemicalEngineMaps::species_upper_bounds)
            .def_property_readonly("species_lower_bounds", &ChemicalEngineMaps::species_lower_bounds)
            .def_property_readonly("phase_species_ln_activities", &ChemicalEngineMaps::phase_species_ln_activities)
            .def_property_readonly("phase_species_ln_activity_coefficients", &ChemicalEngineMaps::phase_species_ln_activity_coefficients)
            .def_property_readonly("phase_species_upper_bounds", &ChemicalEngineMaps::phase_species_upper_bounds)
            .def_property_readonly("phase_species_lower_bounds", &ChemicalEngineMaps::phase_species_lower_bounds)
            .def_property_readonly("solids_mass_frac", &ChemicalEngineMaps::solids_mass_frac)
            .def_property_readonly("solids_volume_frac", &ChemicalEngineMaps::solids_volume_frac)
            .def_property_readonly("aq_volume_frac", &ChemicalEngineMaps::aq_volume_frac)
            .def_property_readonly("phases_volume", &ChemicalEngineMaps::phases_volume)
            .def_property_readonly("phases_mass", &ChemicalEngineMaps::phases_mass)
            .def_property_readonly("phases_volume_frac", &ChemicalEngineMaps::phases_volume_frac)

            .def("add_multiple_species_amt", &ChemicalEngineMaps::add_multiple_species_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_species_amt", &ChemicalEngineMaps::add_species_amt, py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("add_element_amt", &ChemicalEngineMaps::add_element_amt, py::arg("element_name"), py::arg("val"), py::arg("units")="moles")
            .def("add_multiple_elements_amt", &ChemicalEngineMaps::add_multiple_elements_amt, py::arg("input_dict"), py::arg("units")="moles")
            .def("add_amt_from_formula", &ChemicalEngineMaps::add_amt_from_formula, py::arg("formula"), py::arg("val"), py::arg("units")="moles")
            .def("get_b_from_formula", &ChemicalEngineMaps::get_b_from_formula, py::arg("formula"), py::arg("val")=1, py::arg("units")="moles", py::arg("min_amount")=1e-15)

            .def("set_multiple_species_lower_bound", &ChemicalEngineMaps::set_multiple_species_lower_bound, py::arg("input_dict"), py::arg("units")="moles")
            .def("set_multiple_species_upper_bound", &ChemicalEngineMaps::set_multiple_species_upper_bound, py::arg("input_dict"), py::arg("units")="moles")
            .def("set_species_lower_bound", py::overload_cast<const std::string&, double, const std::string&>
                 (&ChemicalEngineMaps::set_species_lower_bound), py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_lower_bound", py::overload_cast<Index, double, const std::string&>
                 (&ChemicalEngineMaps::set_species_lower_bound), py::arg("ispecies"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_upper_bound", py::overload_cast<const std::string&, double, const std::string&>
                 (&ChemicalEngineMaps::set_species_upper_bound), py::arg("species"), py::arg("val"), py::arg("units")="moles")
            .def("set_species_upper_bound", py::overload_cast<Index, double, const std::string&>
                 (&ChemicalEngineMaps::set_species_upper_bound), py::arg("ispecies"), py::arg("val"), py::arg("units")="moles")

            .def("supress_phase", &ChemicalEngineMaps::supress_phase)
            .def("supress_multiple_phases", &ChemicalEngineMaps::supress_multiple_phases)
            .def("supress_species", &ChemicalEngineMaps::supress_species)
            .def("supress_multiple_species", &ChemicalEngineMaps::supress_multiple_species)
            .def("activate_phase", &ChemicalEngineMaps::activate_phase)
            .def("activate_multiple_phases", &ChemicalEngineMaps::activate_multiple_phases)
            .def("activate_species", &ChemicalEngineMaps::activate_species)
            .def("activate_multiple_species", &ChemicalEngineMaps::activate_multiple_species)
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
