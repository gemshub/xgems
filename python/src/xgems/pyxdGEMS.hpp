// ChemicalFun is a C++ and Python library for of C++ and Python API
// for Chemical Formula Parser and Reactions Generator.
//
// Copyright (C) 2018-2026 G.D.Miron, D.Kulik, S.Dmytriieva
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
Gems interface in calculator format for easy using dictionaries.
)doc");

    gems.def(py::init<std::string, bool, bool>(),
             py::arg("input_file"), py::arg("reset_calc")=false, py::arg("cold_start")=true,
             R"doc(
Constructs a ChemicalEngineDicts instance by loading a GEM-Selektor project file.

:param str input_file: The file path for the chemical system definition (e.g., "my-system-dat.lst").
:param bool reset_calc: If true, clear the amounts of all elements, default false.
:param bool cold_start: If true, configures the engine to use a cold start, default true.

**Example:**

.. code-block:: python

    engine = ChemicalEngineDicts("my-system-dat.lst")
)doc" )
            .def(py::init<std::string, std::string, std::string, bool, bool>(),
             py::arg("dch_json"), py::arg("ipm_json"), py::arg("dbr_json"),
             py::arg("reset_calc")=false, py::arg("cold_start")=true,
             R"doc(
Constructs a ChemicalEngineDicts instance from three JSON strings.

:param str dch_json: JSON string for the chemical system definition (*-dch.json content).
:param str ipm_json: JSON string for IPM parameters (*-ipm.json content).
:param str dbr_json: JSON string for node bulk composition (*-dbr.json content).
:param bool reset_calc: If true, clear the amounts of all elements, default false.
:param bool cold_start: If true, use cold start, default true.

**Example:**

.. code-block:: python

    engine = ChemicalEngineDicts(dch_json, ipm_json, dbr_json)
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
Property pressure in Pascals (Pa).

**Example:**

.. code-block:: python

    engine.P = 100000.0
)doc")
            .def_property_readonly("element_names", &ChemicalEngineMaps::element_names,
             R"doc(
Read-only property: the names of all elements in the system.

:return list[str]: List of elements in the system.

**Example:**

.. code-block:: python

    print("element_names", engine.element_names)
)doc")
            .def_property_readonly("species_names", &ChemicalEngineMaps::species_names,
             R"doc(
Read-only property: the names of all species in the system.

:return list[str]: List of species in the system.

**Example:**

.. code-block:: python

    print("species_names", engine.species_names)
)doc")
            .def_property_readonly("phase_names", &ChemicalEngineMaps::phase_names,
             R"doc(
Read-only property: the names of all phases in the system.

:return list[str]: List of phases in the system.

**Example:**

.. code-block:: python

    print("phase_names", engine.phase_names)
)doc")
            .def_property_readonly("aq_phase_symbol", &ChemicalEngineMaps::aq_phase_symbol,
             R"doc(
Read-only property: the aqueous phase name.

:return str: Aqueous phase name. If empty, the aqueous phase is not in system.

**Example:**

.. code-block:: python

    print("aq_phase_symbol", engine.aq_phase_symbol)
)doc")
            .def_property_readonly("gas_phase_symbol", &ChemicalEngineMaps::gas_phase_symbol,
             R"doc(
Read-only property: the gaseous phase name.

:return str: Gaseous phase name. If empty, the gaseous phase is not in system.

**Example:**

.. code-block:: python

    print("gas_phase_symbol", engine.gas_phase_symbol)
)doc")
            .def_property_readonly("gas_species_partial_pressures", &ChemicalEngineMaps::gas_species_partial_pressures,
             R"doc(
Read-only property: partial pressures of all gas species in Pa (x_i · P_total).

Returns an empty dict if no gas phase is present.

:return dict[str, float]: Gas species partial pressures in Pa.

**Example:**

.. code-block:: python

    pp = engine.gas_species_partial_pressures
    print("CO2 partial pressure:", pp.get("CO2"), "Pa")
)doc")
            .def_property_readonly("gas_species_fugacities", &ChemicalEngineMaps::gas_species_fugacities,
             R"doc(
Read-only property: fugacities of all gas species in Pa (exp(ln a_i) · P°, P° = 1 bar).

Returns an empty dict if no gas phase is present.

:return dict[str, float]: Gas species fugacities in Pa.

**Example:**

.. code-block:: python

    f = engine.gas_species_fugacities
    print("CO2 fugacity:", f.get("CO2"), "Pa")
)doc")
            .def_property_readonly("gas_species_fugacity_coefficients", &ChemicalEngineMaps::gas_species_fugacity_coefficients,
             R"doc(
Read-only property: fugacity coefficients of all gas species (dimensionless, exp(ln γ_i)).

Returns an empty dict if no gas phase is present.

:return dict[str, float]: Gas species fugacity coefficients.

**Example:**

.. code-block:: python

    phi = engine.gas_species_fugacity_coefficients
    print("CO2 fugacity coefficient:", phi.get("CO2"))
)doc")
            .def("gas_partial_pressure", &ChemicalEngineMaps::gas_partial_pressure,
             R"doc(
Partial pressure of a single gas species (Pa).

:param str species: Species name in the gas phase.
:return float: Partial pressure in Pa, or 0 if not in gas phase.

**Example:**

.. code-block:: python

    pCO2 = engine.gas_partial_pressure("CO2")
)doc")
            .def("gas_fugacity", &ChemicalEngineMaps::gas_fugacity,
             R"doc(
Fugacity of a single gas species (Pa).

:param str species: Species name in the gas phase.
:return float: Fugacity in Pa, or 0 if not in gas phase.

**Example:**

.. code-block:: python

    fCO2 = engine.gas_fugacity("CO2")
)doc")
            .def("gas_fugacity_coefficient", &ChemicalEngineMaps::gas_fugacity_coefficient,
             R"doc(
Fugacity coefficient of a single gas species (dimensionless).

:param str species: Species name in the gas phase.
:return float: Fugacity coefficient, or 0 if not in gas phase.

**Example:**

.. code-block:: python

    phi_CO2 = engine.gas_fugacity_coefficient("CO2")
)doc")
            .def_property_readonly("element_molar_masses", &ChemicalEngineMaps::element_molar_masses,
             R"doc(
Read-only property: the dictionary of molar masses (kg/mol) for each element.

:return dict[str:float]: Molar masses (kg/mol) for each element.

**Example:**

.. code-block:: python

    print("element_molar_masses", engine.element_molar_masses)
)doc")
            .def_property_readonly("species_in_phase", &ChemicalEngineMaps::species_in_phase,
             R"doc(
Read-only property: the dictionary of the names of all species for each phase in the system.

:return dict[str:list[str]]: List of all species for each phase in the system.

**Example:**

.. code-block:: python

    out = engine.species_in_phase
    for name in engine.phase_names:  #print fixed order to compare
        print(name, out[name])
)doc")
            .def_property_readonly("species_charges", &ChemicalEngineMaps::species_charges,
             R"doc(
Read-only property: the dictionary of the electrical charge of a species.

:return dict[str:float]: The electrical charge of a species.

**Example:**

.. code-block:: python

    out = engine.species_charges
    for name in engine.species_names:  #print fixed order to compare
        print(name, out[name])
)doc")
            .def_property_readonly("species_molar_mass", &ChemicalEngineMaps::species_molar_mass,
             R"doc(
Read-only property: the dictionary of the species molar masses (kg/mol).

:return dict[str:float]: The species molar masses (kg/mol).

**Example:**

.. code-block:: python

    out = engine.species_molar_mass
    for name in engine.species_names:  #print fixed order to compare
        print(name, out[name])
)doc")
            .def_property_readonly("species_molar_volumes", &ChemicalEngineMaps::species_molar_volumes,
             R"doc(
Read-only property: the dictionary of the species standard molar volumes in m³/mol.

:return dict[str:float]: The species standard molar volumes in m³/mol.

**Example:**

.. code-block:: python

    out = engine.species_molar_volumes
    for name in engine.species_names:  #print fixed order to compare
        print(name, out[name])
)doc")

            .def("equilibrate", static_cast<std::string(ChemicalEngineMaps::*)()>(&ChemicalEngineMaps::equilibrate),
             R"doc(
Computes the equilibrium state of the current system.
Uses current temperature (K), pressure (Pa), and element amounts (in mol) to compute equilibrium.

:return str: The string indicating the status.

**Example:**

.. code-block:: python

    engine.T = 298.15
    engine.P = 100000.0
    bulk_composition = {'C': 1e-08, 'Ca': 1e-08, 'Cl': 0.002, 'H': 111.016746657646,
                        'Mg': 0.001, 'O': 55.5083933588231, 'Sn': 130.841288437146, 'Zz': 0.0}
    engine.set_bulk_composition(bulk_composition)
    engine.equilibrate()
)doc")
            .def("equilibrate",
                 static_cast<std::string(ChemicalEngineMaps::*)(double, double, ValuesMap, double)>(&ChemicalEngineMaps::equilibrate),
                 py::arg("T"), py::arg("P"), py::arg("b_dict"), py::arg("min_amount")=1e-15,
             R"doc(
Computes the equilibrium state with explicit temperature, pressure, and bulk composition.

Mirrors ``ChemicalEngine.equilibrate(T, P, b)`` but accepts a dictionary for the bulk
composition. Sets the internal T, P, and bulk composition before computing equilibrium.

:param float T: Temperature in Kelvin.
:param float P: Pressure in Pascals.
:param dict b_dict: Dictionary of element amounts in mol (elements not listed keep their current value).
:param float min_amount: Minimum amount for unspecified elements, default 1e-15.
:return str: The string indicating the status.

**Example:**

.. code-block:: python

    b = {'C': 0.001, 'Ca': 1e-9, 'Cl': 0.002, 'H': 110.68, 'Mg': 0.001, 'O': 55.34, 'Zz': 0.0}
    result = engine.equilibrate(298.15, 101325.0, b)
    print(result)
)doc")
            .def("reequilibrate", static_cast<std::string(ChemicalEngineMaps::*)()>(&ChemicalEngineMaps::reequilibrate),
             R"doc(
Re-equilibrates the system using the current internal state without changing T, P, or bulk composition.

Mirrors ``ChemicalEngine.reequilibrate()``.

:return str: The string indicating the status.

**Example:**

.. code-block:: python

    result = engine.reequilibrate()
    print(result)
)doc")
            .def("reequilibrate", static_cast<std::string(ChemicalEngineMaps::*)(bool)>(&ChemicalEngineMaps::reequilibrate),
                 py::arg("warmstart"),
             R"doc(
Re-equilibrates the system with explicit warm/cold start control.

Mirrors ``ChemicalEngine.reequilibrate(warmstart)``.

:param bool warmstart: If true, uses previous speciation as initial guess (SIA).
:return str: The string indicating the status.

**Example:**

.. code-block:: python

    result = engine.reequilibrate(True)
    print(result)
)doc")
            .def("setPT", &ChemicalEngineMaps::setPT,
                 py::arg("P"), py::arg("T"),
             R"doc(
Sets the pressure and temperature without computing equilibrium.

Mirrors ``ChemicalEngine.setPT(P, T)``. Also updates the T and P member variables.

:param float P: Pressure in Pascals.
:param float T: Temperature in Kelvin.
:return bool: True if PT was set correctly, False if out of range.

**Example:**

.. code-block:: python

    ok = engine.setPT(101325.0, 298.15)
)doc")
            .def("setB", &ChemicalEngineMaps::setB,
                 py::arg("b_dict"), py::arg("min_amount")=1e-15,
             R"doc(
Sets the amounts of elements from a dictionary.

Mirrors ``ChemicalEngine.setB(b)`` but accepts a dictionary instead of a vector.
Elements absent from ``b_dict`` are left at their current value (or set to ``min_amount``
if they fall below it).

:param dict b_dict: Dictionary of element amounts in mol.
:param float min_amount: Minimum amount for unspecified elements, default 1e-15.

**Example:**

.. code-block:: python

    engine.setB({'Ca': 0.001, 'Cl': 0.002, 'H': 110.0, 'O': 55.0, 'Zz': 0.0})
)doc")
            .def_property_readonly("temperature", &ChemicalEngineMaps::temperature,
             R"doc(
Read-only property: current temperature of the system in K after the last equilibration.

Mirrors ``ChemicalEngine.temperature()``. Use the writable ``T`` attribute to set the
temperature for the next equilibration.

:return float: Temperature in Kelvin.

**Example:**

.. code-block:: python

    print("temperature", engine.temperature)
)doc")
            .def_property_readonly("pressure", &ChemicalEngineMaps::pressure,
             R"doc(
Read-only property: current pressure of the system in Pa after the last equilibration.

Mirrors ``ChemicalEngine.pressure()``. Use the writable ``P`` attribute to set the
pressure for the next equilibration.

:return float: Pressure in Pascals.

**Example:**

.. code-block:: python

    print("pressure", engine.pressure)
)doc")
            .def("readDbrFromFile", &ChemicalEngineMaps::readDbrFromFile,
                 py::arg("filename"),
             R"doc(
Reads a DBR file from disk, updating the system composition.

Mirrors ``ChemicalEngine.readDbrFromFile(filename)``.

:param str filename: Path to the DBR file (e.g., "*-dbr.json").

**Example:**

.. code-block:: python

    engine.readDbrFromFile("system-dbr-0-0000.json")
)doc")
            .def("readDbrFromJsonString", &ChemicalEngineMaps::readDbrFromJsonString,
                 py::arg("dbr_json"),
             R"doc(
Reads the system DBR composition from a JSON string.

Mirrors ``ChemicalEngine.readDbrFromJsonString(dbr_json)``.

:param str dbr_json: JSON string containing DBR composition data.

**Example:**

.. code-block:: python

    engine.readDbrFromJsonString(dbr_json_str)
)doc")
            .def("writeDbrToFile", &ChemicalEngineMaps::writeDbrToFile,
                 py::arg("filename"),
             R"doc(
Writes the current DBR to a file.

Mirrors ``ChemicalEngine.writeDbrToFile(filename)``.

:param str filename: Path to the output DBR file.

**Example:**

.. code-block:: python

    engine.writeDbrToFile("system-dbr-0-0001.dat")
)doc")
            .def("writeDbrToJsonString", &ChemicalEngineMaps::writeDbrToJsonString,
             R"doc(
Returns the current DBR as a JSON string.

Mirrors ``ChemicalEngine.writeDbrToJsonString()``.

:return str: JSON representation of the current system state.

**Example:**

.. code-block:: python

    json_str = engine.writeDbrToJsonString()
    print(json_str)
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

:return int: Number of elements in the system.

**Example:**

.. code-block:: python

    print("nelements", engine.nelements)
)doc")
            .def_property_readonly("nphases", &ChemicalEngineMaps::nphases,
             R"doc(
Read-only property: the number of phases in the system.

:return int: Number of phases in the system.

**Example:**

.. code-block:: python

    print("nphases", engine.nphases)
)doc")
            .def_property_readonly("nspecies", &ChemicalEngineMaps::nspecies,
             R"doc(
Read-only property: the number of species in the system.

:return int: Number of species in the system.

**Example:**

.. code-block:: python

    print("nspecies", engine.nspecies)
)doc")
            .def("clear", &ChemicalEngineMaps::clear, py::arg("min_amount")=1e-15,
             R"doc(
Clear the amounts of elements (set the default amount for all components).

:param float min_amount: Min amount of element in mol, default 1e-15.

**Example:**

.. code-block:: python

    engine.clear(1e-12)
)doc")

                .def("set_species_G0", &ChemicalEngineMaps::set_species_G0,
             py::arg("name"), py::arg("value"), py::arg("phase") = py::none(),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).

:param str name: Species name.
:param float value: Standard molar Gibbs energy value (J/mol).
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

    engine.set_species_G0("H2O", -237.13)
)doc")
            .def_property_readonly("bulk_composition", &ChemicalEngineMaps::bulk_composition,
             R"doc(
Read-only property: the dictionary of the amounts of the elements in mol.

:return dict[str:float]: The amounts of the elements in mol.

**Example:**

.. code-block:: python

    print("bulk_composition", engine.bulk_composition)
)doc")
            .def_property_readonly("pH", &ChemicalEngineMaps::pH,
             R"doc(
Read-only property: the pH of the aqueous phase (in the activity scale (-log10 molal)).

:return float: pH of the aqueous phase.

**Example:**

.. code-block:: python

    print("pH", engine.pH)
)doc")
            .def_property_readonly("pE", &ChemicalEngineMaps::pE,
             R"doc(
Read-only property: the pe of the aqueous phase (in the activity scale (-log10 molal)).

:return float: pe of the aqueous phase.

**Example:**

.. code-block:: python

    print("pE", engine.pE)
)doc")
            .def_property_readonly("ionic_strength", &ChemicalEngineMaps::ionic_strength,
             R"doc(
Read-only property: the ionic strength of the aqueous phase in molal.

:return float: Ionic strength of the aqueous phase.

**Example:**

.. code-block:: python

    print("ionic_strength", engine.ionic_strength)
)doc")
            .def_property_readonly("system_volume", &ChemicalEngineMaps::system_volume,
             R"doc(
Read-only property: the total volume of the system in m³.

:return float: Total volume of the system in m³.

**Example:**

.. code-block:: python

    print("system_volume", engine.system_volume)
)doc")
            .def_property_readonly("system_mass", &ChemicalEngineMaps::system_mass,
             R"doc(
Read-only property: the total mass of the system in kg.

:return float: Total mass of the system in kg.

**Example:**

.. code-block:: python

    print("system_mass", engine.system_mass)
)doc")
            .def_property_readonly("Eh", &ChemicalEngineMaps::Eh,
             R"doc(
Read-only property: the Eh (redox potential) of the aqueous phase in V.

:return float: Eh of the aqueous phase in Volts.

**Example:**

.. code-block:: python

    print("Eh", engine.Eh)
)doc")
            .def_property_readonly("water_activity", &ChemicalEngineMaps::water_activity,
             R"doc(
Read-only property: activity of water (H2O@) in the aqueous phase.

Returns 1.0 if no aqueous phase is present.

:return float: Water activity (dimensionless).

**Example:**

.. code-block:: python

    print("water activity:", engine.water_activity)
)doc")
            .def_property_readonly("osmotic_coefficient", &ChemicalEngineMaps::osmotic_coefficient,
             R"doc(
Read-only property: osmotic coefficient Φ of the aqueous phase.

Computed as Φ = −ln(a_w) × n_water / n_solutes.
Returns 1.0 for pure water or if no aqueous phase is present.

:return float: Osmotic coefficient (dimensionless).

**Example:**

.. code-block:: python

    print("osmotic coefficient:", engine.osmotic_coefficient)
)doc")
            .def_property_readonly("hardness", &ChemicalEngineMaps::hardness,
             R"doc(
Read-only property: total, Ca, and Mg hardness of the aqueous phase (mg/L as CaCO₃).

Uses total Ca and Mg element amounts in the aqueous phase (all species).
Dictionary keys: "total", "Ca", "Mg".

:return dict[str, float]: Hardness components in mg/L as CaCO₃.

**Example:**

.. code-block:: python

    h = engine.hardness
    print("Total hardness:", h["total"], "mg/L as CaCO3")
    print("Ca hardness:", h["Ca"], "mg/L as CaCO3")
)doc")
            .def_property_readonly("converged", &ChemicalEngineMaps::converged,
             R"doc(
Read-only property: whether the last equilibrium computation converged.

:return bool: True if the last equilibration converged.

**Example:**

.. code-block:: python

    print("converged", engine.converged)
)doc")
            .def_property_readonly("mass_balance_errors", &ChemicalEngineMaps::mass_balance_errors,
             R"doc(
Read-only property: per-element mass balance residual (mol) after equilibration.

Computes W·n − b for each element, where W is the formula matrix, n is the
equilibrium species amounts, and b is the input bulk composition. Values near
zero confirm a well-converged result; large values indicate solver failure.

:return dict: Dictionary mapping element names to residuals in mol.

**Example:**

.. code-block:: python

    errors = engine.mass_balance_errors
    print("Residual for Ca:", errors["Ca"], "mol")
)doc")
            .def_property_readonly("mass_balance_relative_errors", &ChemicalEngineMaps::mass_balance_relative_errors,
             R"doc(
Read-only property: per-element relative mass balance residual after equilibration.

Computes (W·n − b) / b for each element. Elements with b ≈ 0 (e.g., charge row
Zz) are returned as 0.

:return dict: Dictionary mapping element names to dimensionless relative residuals.

**Example:**

.. code-block:: python

    rel = engine.mass_balance_relative_errors
    print("Relative residual for Ca:", rel["Ca"])
)doc")
            .def("phase_is_present", &ChemicalEngineMaps::phase_is_present,
             py::arg("phase_name"), py::arg("threshold") = 1e-12,
             R"doc(
Return True if a phase has a molar amount above the threshold.

:param str phase_name: Phase name.
:param float threshold: Minimum amount in mol (default 1e-12).
:return bool: True if the phase is present.

**Example:**

.. code-block:: python

    if engine.phase_is_present("Calcite"):
        print("Calcite is precipitating")
)doc")
            .def("present_phases", &ChemicalEngineMaps::present_phases,
             py::arg("threshold") = 1e-12,
             R"doc(
Return names of all phases whose molar amount exceeds the threshold.

:param float threshold: Minimum amount in mol (default 1e-12).
:return list[str]: Names of present phases.

**Example:**

.. code-block:: python

    for name in engine.present_phases():
        print(name)
)doc")
            .def("present_minerals", &ChemicalEngineMaps::present_minerals,
             py::arg("threshold") = 1e-12,
             R"doc(
Return names of all mineral (non-aqueous, non-gas) phases above the threshold.

:param float threshold: Minimum amount in mol (default 1e-12).
:return list[str]: Names of present mineral phases.

**Example:**

.. code-block:: python

    minerals = engine.present_minerals()
    print("Precipitated minerals:", minerals)
)doc")
            .def("is_aqueous_phase", &ChemicalEngineMaps::is_aqueous_phase,
             R"doc(
Return True if the given phase name is the aqueous phase.

:param str phase: Phase name.
:return bool: True if aqueous.

**Example:**

.. code-block:: python

    engine.is_aqueous_phase("aq_gen")  # True
)doc")
            .def("is_gas_phase", &ChemicalEngineMaps::is_gas_phase,
             R"doc(
Return True if the given phase name is the gas phase.

:param str phase: Phase name.
:return bool: True if gas.

**Example:**

.. code-block:: python

    engine.is_gas_phase("gas_gen")  # True
)doc")
            .def("is_mineral_phase", &ChemicalEngineMaps::is_mineral_phase,
             R"doc(
Return True if the given phase is a mineral (not aqueous, not gas).

:param str phase: Phase name.
:return bool: True if mineral.

**Example:**

.. code-block:: python

    engine.is_mineral_phase("Calcite")  # True
)doc")
            .def("log_K", &ChemicalEngineMaps::log_K,
             R"doc(
Standard log10 equilibrium constant at current T and P.

The reaction is given as a dict mapping species names to stoichiometric
coefficients (positive = product, negative = reactant).

:param dict reaction: {species_name: stoich_coeff}.
:return float: log10(K).

**Example:**

.. code-block:: python

    logK = engine.log_K({"Ca+2": 1, "CO3-2": 1, "Calcite": -1})
)doc")
            .def("delta_G0_reaction", &ChemicalEngineMaps::delta_G0_reaction,
             R"doc(
Standard molar Gibbs energy change of a reaction at current T and P (J/mol).

:param dict reaction: {species_name: stoich_coeff}.
:return float: ΔG° in J/mol.

**Example:**

.. code-block:: python

    dG = engine.delta_G0_reaction({"Ca+2": 1, "CO3-2": 1, "Calcite": -1})
)doc")
            .def("delta_H0_reaction", &ChemicalEngineMaps::delta_H0_reaction,
             R"doc(
Standard molar enthalpy change of a reaction at current T and P (J/mol).

:param dict reaction: {species_name: stoich_coeff}.
:return float: ΔH° in J/mol.

**Example:**

.. code-block:: python

    dH = engine.delta_H0_reaction({"Ca+2": 1, "CO3-2": 1, "Calcite": -1})
)doc")
            .def_property_readonly("num_iterations", &ChemicalEngineMaps::num_iterations,
             R"doc(
Read-only property: number of iterations of the last equilibrium computation.

:return int: Number of iterations.

**Example:**

.. code-block:: python

    print("num_iterations", engine.num_iterations)
)doc")
            .def_property_readonly("elapsed_time", &ChemicalEngineMaps::elapsed_time,
             R"doc(
Read-only property: elapsed time of the last equilibrium computation in seconds.

:return float: Elapsed time in seconds.

**Example:**

.. code-block:: python

    print("elapsed_time", engine.elapsed_time)
)doc")
            .def_property_readonly("system_gibbs_energy", &ChemicalEngineMaps::system_gibbs_energy,
             R"doc(
Read-only property: the total Gibbs energy of the system in J.

:return float: Total Gibbs energy in J.

**Example:**

.. code-block:: python

    print("system_gibbs_energy", engine.system_gibbs_energy)
)doc")
            .def_property_readonly("system_enthalpy", &ChemicalEngineMaps::system_enthalpy,
             R"doc(
Read-only property: the total enthalpy of the system in J.

:return float: Total enthalpy in J.

**Example:**

.. code-block:: python

    print("system_enthalpy", engine.system_enthalpy)
)doc")
            .def_property_readonly("system_entropy", &ChemicalEngineMaps::system_entropy,
             R"doc(
Read-only property: the total entropy of the system in J/K.

:return float: Total entropy in J/K.

**Example:**

.. code-block:: python

    print("system_entropy", engine.system_entropy)
)doc")
            .def_property_readonly("system_heat_capacity_const_p", &ChemicalEngineMaps::system_heat_capacity_const_p,
             R"doc(
Read-only property: the total isobaric heat capacity of the system in J/K.

:return float: Total isobaric heat capacity in J/K.

**Example:**

.. code-block:: python

    print("system_heat_capacity_const_p", engine.system_heat_capacity_const_p)
)doc")
            .def_property_readonly("phases_molar_volume", &ChemicalEngineMaps::phases_molar_volume,
             R"doc(
Read-only property: the dictionary of the phases molar volumes in m³/mol.

:return dict[str:float]: The phases molar volumes in m³/mol.

**Example:**

.. code-block:: python

    print("phases_molar_volume", engine.phases_molar_volume)
)doc")
            .def_property_readonly("phase_sat_indices", &ChemicalEngineMaps::phase_sat_indices,
             R"doc(
Read-only property: the dictionary of the saturation indices of all phases (log₁₀ units).

:return dict[str:float]: The saturation indices of all phases.

**Example:**

.. code-block:: python

    print("phase_sat_indices", engine.phase_sat_indices)
)doc")
            .def_property_readonly("phases_density", &ChemicalEngineMaps::phases_density,
             R"doc(
Read-only property: the dictionary of the densities of all phases in kg/m³.

:return dict[str:float]: The densities of all phases in kg/m³.

**Example:**

.. code-block:: python

    print("phases_density", engine.phases_density)
)doc")
            .def_property_readonly("phases_enthalpy", &ChemicalEngineMaps::phases_enthalpy,
             R"doc(
Read-only property: the dictionary of the total enthalpies of all phases in J.

:return dict[str:float]: The total enthalpies of all phases in J.

**Example:**

.. code-block:: python

    print("phases_enthalpy", engine.phases_enthalpy)
)doc")
            .def_property_readonly("phases_entropy", &ChemicalEngineMaps::phases_entropy,
             R"doc(
Read-only property: the dictionary of the total entropies of all phases in J/K.

:return dict[str:float]: The total entropies of all phases in J/K.

**Example:**

.. code-block:: python

    print("phases_entropy", engine.phases_entropy)
)doc")
            .def_property_readonly("phases_molar_gibbs_energy", &ChemicalEngineMaps::phases_molar_gibbs_energy,
             R"doc(
Read-only property: the dictionary of the molar Gibbs energies of all phases in J/mol.

:return dict[str:float]: The molar Gibbs energies of all phases in J/mol.

**Example:**

.. code-block:: python

    print("phases_molar_gibbs_energy", engine.phases_molar_gibbs_energy)
)doc")
            .def_property_readonly("phases_heat_capacity_const_p", &ChemicalEngineMaps::phases_heat_capacity_const_p,
             R"doc(
Read-only property: the dictionary of the total isobaric heat capacities of all phases in J/K.

:return dict[str:float]: The total isobaric heat capacities of all phases in J/K.

**Example:**

.. code-block:: python

    print("phases_heat_capacity_const_p", engine.phases_heat_capacity_const_p)
)doc")

            .def_property_readonly("aq_elements_molarity", &ChemicalEngineMaps::aq_elements_molarity,
             R"doc(
Read-only property: the dictionary for aq elements the aq solution composition in mol/L aq solution.

:return dict[str:float]: The aq solution composition in mol/L.

**Example:**

.. code-block:: python

    print("aq_elements_molarity", engine.aq_elements_molarity)
)doc")
            .def_property_readonly("aq_elements_molality", &ChemicalEngineMaps::aq_elements_molality,
             R"doc(
Read-only property: the dictionary for aq elements the aq solution elemental composition in mol/kgH2O.

:return dict[str:float]: The aq solution elemental composition in mol/kgH2O.

**Example:**

.. code-block:: python

    print("aq_elements_molality", engine.aq_elements_molality)
)doc")
            .def_property_readonly("aq_species_molarity", &ChemicalEngineMaps::aq_species_molarity,
             R"doc(
Read-only property: the dictionary for aq species the aq solution composition in mol/L of aqueous species.

:return dict[str:float]: The aq solution composition in mol/L of aqueous species.

**Example:**

.. code-block:: python

    print("aq_species_molarity", engine.aq_species_molarity)
)doc")
            .def_property_readonly("aq_species_molality", &ChemicalEngineMaps::aq_species_molality,
             R"doc(
Read-only property: the dictionary for the aq solution composition in mol/kg H2O of aqueous species (speciation).

:return dict[str:float]: The aq solution composition in mol/kg H2O of aqueous species.

**Example:**

.. code-block:: python

    print("aq_species_molality", engine.aq_species_molality)
)doc")
            .def_property_readonly("aq_elements_moles", &ChemicalEngineMaps::aq_elements_moles,
             R"doc(
Read-only property: the dictionary of the amounts of each element in the aqueous phase (in mol).

:return dict[str:float]: The amounts of each element in the aqueous phase (in mol).

**Example:**

.. code-block:: python

    print("aq_elements_moles", engine.aq_elements_moles)
)doc")
            .def("set_bulk_composition", &ChemicalEngineMaps::set_bulk_composition,
             py::arg("b_input"), py::arg("min_amount")=1e-15,
             R"doc(
Sets the amounts of elements (vector b).

:param dict b_input: Dictionary of elements amounts in mol.
:param float min_amount: Min amount of element in mol, default 1e-15.

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

:param float min_amount: Min amount of element in mol, default 1e-15.

**Example:**

.. code-block:: python

    engine.reset_aq_composition()
)doc")
            .def("solids_elements_moles", &ChemicalEngineMaps::solids_elements_moles,
             py::arg("min_amount_phase")=1e-12, py::arg("min_amount_element")=1e-14,
             R"doc(
Returns the dictionary of the mole amounts of elements in all solids together.

:param float min_amount_phase: Min amount of phase in mol, default 1e-12.
:param float min_amount_element: Min amount of element in mol, default 1e-14.
:return dict[str:float]: The mole amounts of elements.

**Example:**

.. code-block:: python

    print("solids_elements_moles", engine.solids_elements_moles())
)doc")
            .def_property_readonly("phases_elements_moles", &ChemicalEngineMaps::phases_elements_moles,
             R"doc(
Read-only property: the dictionary of dictionaries containing mole amounts of elements for each phase (in mol).

:return dict[str:dict[str:float]]: The mole amounts of elements for each phase (in mol).

**Example:**

.. code-block:: python

    print("phases_elements_moles", engine.phases_elements_moles)
)doc")

            .def("phase_species_moles", py::overload_cast<>(&ChemicalEngineMaps::phase_species_moles),
             R"doc(
Get the dictionary of dictionaries containing species amounts in mol for each phase.

:return dict[str:dict[str:float]]: The species amounts in mol for each phase.

**Example:**

.. code-block:: python

    out = engine.phase_species_moles()
    for name in engine.phase_names:
        print(name, out[name])
)doc")
            .def("phase_species_moles", py::overload_cast<std::string>(&ChemicalEngineMaps::phase_species_moles),
             R"doc(
Get the dictionary of phase species amounts in mol.

:param str name: Phase name.
:return dict[str:float]: The species amounts in mol.

**Example:**

.. code-block:: python

    print("phase_species_moles ('aq_gen')", engine.phase_species_moles('aq_gen'))
)doc")
            .def_property_readonly("phases_moles", &ChemicalEngineMaps::phases_moles,
             R"doc(
Read-only property: the dictionary of the molar amounts of all phases in mol.

:return dict[str:float]: The molar amounts of all phases in mol.

**Example:**

.. code-block:: python

    print("phases_moles", engine.phases_moles)
)doc")
            .def_property_readonly("species_moles", &ChemicalEngineMaps::species_moles,
             R"doc(
Read-only property: the dictionary of the amounts of all species in mol.

:return dict[str:float]: The amounts of all species in mol.

**Example:**

.. code-block:: python

    print("species_moles", engine.species_moles)
)doc")
            .def_property_readonly("species_ln_activities", &ChemicalEngineMaps::species_ln_activities,
             R"doc(
Read-only property: the dictionary of the ln activities of all species.

:return dict[str:float]: The ln activities of all species.

**Example:**

.. code-block:: python

    print("species_ln_activities", engine.species_ln_activities)
)doc")
            .def_property_readonly("species_ln_activity_coefficients", &ChemicalEngineMaps::species_ln_activity_coefficients,
             R"doc(
Read-only property: the dictionary of the ln activity coefficients of all species (mole fraction scale).

:return dict[str:float]: The ln activity coefficients of all species.

**Example:**

.. code-block:: python

    print("species_ln_activity_coefficients", engine.species_ln_activity_coefficients)
)doc")
            .def_property_readonly("species_mole_fractions", &ChemicalEngineMaps::species_mole_fractions,
             R"doc(
Read-only property: the dictionary of the mole fractions of all species.

:return dict[str:float]: The mole fractions of all species.

**Example:**

.. code-block:: python

    print("species_mole_fractions", engine.species_mole_fractions)
)doc")
            .def_property_readonly("species_ln_concentrations", &ChemicalEngineMaps::species_ln_concentrations,
             R"doc(
Read-only property: the dictionary of the ln concentrations of all species.

:return dict[str:float]: The ln concentrations of all species.

**Example:**

.. code-block:: python

    print("species_ln_concentrations", engine.species_ln_concentrations)
)doc")
            .def_property_readonly("species_chemical_potentials", &ChemicalEngineMaps::species_chemical_potentials,
             R"doc(
Read-only property: the dictionary of the chemical potentials of all species in J/mol.

:return dict[str:float]: The chemical potentials of all species in J/mol.

**Example:**

.. code-block:: python

    print("species_chemical_potentials", engine.species_chemical_potentials)
)doc")
            .def_property_readonly("species_gibbs_energies", &ChemicalEngineMaps::species_gibbs_energies,
             R"doc(
Read-only property: the dictionary of the standard molar Gibbs energies of all species in J/mol.

:return dict[str:float]: The standard molar Gibbs energies in J/mol.

**Example:**

.. code-block:: python

    print("species_gibbs_energies", engine.species_gibbs_energies)
)doc")
            .def_property_readonly("species_enthalpies", &ChemicalEngineMaps::species_enthalpies,
             R"doc(
Read-only property: the dictionary of the standard molar enthalpies of all species in J/mol.

:return dict[str:float]: The standard molar enthalpies in J/mol.

**Example:**

.. code-block:: python

    print("species_enthalpies", engine.species_enthalpies)
)doc")
            .def_property_readonly("species_entropies", &ChemicalEngineMaps::species_entropies,
             R"doc(
Read-only property: the dictionary of the standard molar entropies of all species in J/mol/K.

:return dict[str:float]: The standard molar entropies in J/mol/K.

**Example:**

.. code-block:: python

    print("species_entropies", engine.species_entropies)
)doc")
            .def_property_readonly("species_heat_capacities_const_p", &ChemicalEngineMaps::species_heat_capacities_const_p,
             R"doc(
Read-only property: the dictionary of the standard molar isobaric heat capacities of all species in J/mol/K.

:return dict[str:float]: The standard molar heat capacities in J/mol/K.

**Example:**

.. code-block:: python

    print("species_heat_capacities_const_p", engine.species_heat_capacities_const_p)
)doc")
            .def_property_readonly("species_upper_bounds", &ChemicalEngineMaps::species_upper_bounds,
             R"doc(
Read-only property: the dictionary of the species upper limits in mol.

:return dict[str:float]: The species upper limits in mol.

**Example:**

.. code-block:: python

    out = engine.species_upper_bounds
    for name in engine.species_names:
        print(name, out[name])
)doc")
            .def_property_readonly("species_lower_bounds", &ChemicalEngineMaps::species_lower_bounds,
             R"doc(
Read-only property: the dictionary of the species lower limits in mol.

:return dict[str:float]: The species lower limits in mol.

**Example:**

.. code-block:: python

    out = engine.species_lower_bounds
    for name in engine.species_names:
        print(name, out[name])
)doc")
            .def_property_readonly("phase_species_ln_activities", &ChemicalEngineMaps::phase_species_ln_activities,
             R"doc(
Read-only property: the dictionary of dictionaries containing species ln activities for each phase.

:return dict[str:dict[str:float]]: The species ln activities for each phase.

**Example:**

.. code-block:: python

    print("phase_species_ln_activities", engine.phase_species_ln_activities)
)doc")
            .def_property_readonly("phase_species_ln_activity_coefficients", &ChemicalEngineMaps::phase_species_ln_activity_coefficients,
             R"doc(
Read-only property: the dictionary of dictionaries containing species ln activity coefficients for each phase.

:return dict[str:dict[str:float]]: The species ln activity coefficients for each phase.

**Example:**

.. code-block:: python

    out = engine.phase_species_ln_activity_coefficients
    for name in engine.phase_names:
        print(name, out[name])
)doc")
            .def_property_readonly("phase_species_upper_bounds", &ChemicalEngineMaps::phase_species_upper_bounds,
             R"doc(
Read-only property: the dictionary of dictionaries containing species upper limits in mol for each phase.

:return dict[str:dict[str:float]]: The species upper limits in mol for each phase.

**Example:**

.. code-block:: python

    print("phase_species_upper_bounds\n", engine.phase_species_upper_bounds)
)doc")
            .def_property_readonly("phase_species_lower_bounds", &ChemicalEngineMaps::phase_species_lower_bounds,
             R"doc(
Read-only property: the dictionary of dictionaries containing species lower limits in mol for each phase.

:return dict[str:dict[str:float]]: The species lower limits in mol for each phase.

**Example:**

.. code-block:: python

    print("phase_species_lower_bounds\n", engine.phase_species_lower_bounds)
)doc")
            .def_property_readonly("solids_mass_frac", &ChemicalEngineMaps::solids_mass_frac,
             R"doc(
Read-only property: the dictionary of the mass(phase)/mass(system) ratios for solid phases.

:return dict[str:float]: The mass(phase)/mass(system) ratios for solid phases.

**Example:**

.. code-block:: python

    print("solids_mass_frac", engine.solids_mass_frac)
)doc")
            .def_property_readonly("solids_volume_frac", &ChemicalEngineMaps::solids_volume_frac,
             R"doc(
Read-only property: the dictionary of the volume(phase)/volume(total) ratio for solid phases.

:return dict[str:float]: The volume(phase)/volume(total) ratio for solid phases.

**Example:**

.. code-block:: python

    print("solids_volume_frac", engine.solids_volume_frac)
)doc")
            .def_property_readonly("aq_volume_frac", &ChemicalEngineMaps::aq_volume_frac,
             R"doc(
Read-only property: the volume fraction of aqueous phase from total system volume.

:return float: The volume fraction of aqueous phase from total system volume.

**Example:**

.. code-block:: python

    print("aq_volume_frac", engine.aq_volume_frac)
)doc")
            .def_property_readonly("phases_volume", &ChemicalEngineMaps::phases_volume,
             R"doc(
Read-only property: the dictionary of phases volumes in m³.

:return dict[str:float]: The phases volumes in m³.

**Example:**

.. code-block:: python

    print("phases_volume", engine.phases_volume)
)doc")
            .def_property_readonly("phases_mass", &ChemicalEngineMaps::phases_mass,
             R"doc(
Read-only property: the dictionary of phases masses in kg.

:return dict[str:float]: The phases masses in kg.

**Example:**

.. code-block:: python

    print("phases_mass", engine.phases_mass)
)doc")
            .def_property_readonly("phases_volume_frac", &ChemicalEngineMaps::phases_volume_frac,
             R"doc(
Read-only property: the dictionary of phases and their volume fractions from total system volume.

:return dict[str:float]: The phases volume fractions from total system volume.

**Example:**

.. code-block:: python

    print("phases_volume_frac", engine.phases_volume_frac)
)doc")

            .def("add_multiple_species_amt", &ChemicalEngineMaps::add_multiple_species_amt,
             py::arg("input_dict"), py::arg("units")="moles",
             R"doc(
Add multiple species amounts in the system useful for adding aqueous solution composition.

:param dict input_dict: Dictionary of species amount in units.
:param str units: Units of amount ("moles", "kg", "m3"), default "moles".

**Example:**

.. code-block:: python

    engine.add_multiple_species_amt({'HCl@': 0.01, 'H2@': 2}, "moles")
)doc")
            .def("add_species_amt", &ChemicalEngineMaps::add_species_amt,
             py::arg("species"), py::arg("val"), py::arg("units")="moles", py::arg("phase") = py::none(),
             R"doc(
Add species amount in the system useful for adding aqueous solution composition.

:param str species: Species symbol.
:param float val: Species amount in units.
:param str units: Units of amount ("moles", "kg", "m3"), default "moles".
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

    engine.add_species_amt('H2O@', 0.01, "kg")
)doc")
            .def("add_element_amt", &ChemicalEngineMaps::add_element_amt,
             py::arg("element_name"), py::arg("val"), py::arg("units")="moles",
             R"doc(
Add element amount in the system.

:param str element_name: Element symbol.
:param float val: Element amount in units.
:param str units: Units of amount ("moles", "kg"), default "moles".

**Example:**

.. code-block:: python

    engine.add_element_amt('Al', 0.3, "moles")
)doc")
            .def("add_multiple_elements_amt", &ChemicalEngineMaps::add_multiple_elements_amt,
             py::arg("input_dict"), py::arg("units")="moles",
             R"doc(
Add multiple elements amount in the system useful for adding aqueous solution composition.

:param dict input_dict: Dictionary of elements amount in units.
:param str units: Units of amount ("moles", "kg"), default "moles".

**Example:**

.. code-block:: python

    engine.add_multiple_elements_amt({'Na': 1.013077, 'Si': 1.013077}, "moles")
)doc")
            .def("add_amt_from_formula", &ChemicalEngineMaps::add_amt_from_formula,
             py::arg("formula"), py::arg("val"), py::arg("units")="moles",
             R"doc(
Add multiple elements using user defined formula.

:param dict formula: User defined formula.
:param float val: Component amount in units.
:param str units: Units of amount ("moles", "kg"), default "moles".

**Example:**

.. code-block:: python

    engine.add_amt_from_formula({'K': 2, 'O': 1}, 4.108e-3, "kg")
)doc")
            .def("get_b_from_formula", &ChemicalEngineMaps::get_b_from_formula,
             py::arg("formula"), py::arg("val")=1, py::arg("units")="moles", py::arg("min_amount")=1e-15,
             R"doc(
Returns a bulk vector b from user-defined formula (as dict {"H":2,"O":1})
and amount of the formula object in units of "moles" or "kg".

:param dict formula: User defined formula.
:param float val: Amount of the formula object in units, default 1.
:param str units: Units of amount ("moles", "kg"), default "moles".
:param float min_amount: Min amount of element in mol, default 1e-15.
:return list[float]: Vector of element amounts in mol.

**Example:**

.. code-block:: python

    b_from_formula = engine.get_b_from_formula({"H": 2, "O": 1}, 1, "kg")
    print("b_from_formula", b_from_formula)
)doc")

            .def("set_multiple_species_lower_bound", &ChemicalEngineMaps::set_multiple_species_lower_bound,
             py::arg("input_dict"), py::arg("units")="moles",
             R"doc(
Sets a lower bound for multiple species.

:param dict input_dict: Dictionary of species lower bound.
:param str units: Units of lower bound ("moles", "kg", "m3"), default "moles".

**Example:**

.. code-block:: python

    engine.set_multiple_species_lower_bound({'Mg(CO3)@': 30, 'Mg(HCO3)+': 40, 'Mg+2': 50})
)doc")
            .def("set_multiple_species_upper_bound", &ChemicalEngineMaps::set_multiple_species_upper_bound,
             py::arg("input_dict"), py::arg("units")="moles",
             R"doc(
Sets an upper bound for multiple species.

:param dict input_dict: Dictionary of species upper bound.
:param str units: Units of upper bound ("moles", "kg", "m3"), default "moles".

**Example:**

.. code-block:: python

    engine.set_multiple_species_upper_bound({'Mg(CO3)@': 300, 'Mg(HCO3)+': 400, 'Mg+2': 500})
)doc")
            .def("set_species_lower_bound", py::overload_cast<const std::string&, double, const std::string&, std::optional<std::string> >
                 (&ChemicalEngineMaps::set_species_lower_bound),
             py::arg("species"), py::arg("val"), py::arg("units")="moles", py::arg("phase") = py::none(),
             R"doc(
Sets a lower bound for a species identified by name.

:param str species: Species name.
:param float val: Lower limit in units.
:param str units: Units of amount ("moles", "kg", "m3"), default "moles".
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

    engine.set_species_lower_bound('Ca(HCO3)+', 200, "moles")
)doc")
            .def("set_species_lower_bound", py::overload_cast<Index, double, const std::string&>
                 (&ChemicalEngineMaps::set_species_lower_bound),
             py::arg("ispecies"), py::arg("val"), py::arg("units")="moles",
             R"doc(
Sets a lower bound for a species identified by its index.

:param int ispecies: Index of the species.
:param float val: Lower limit in units.
:param str units: Units of amount ("moles", "kg", "m3"), default "moles".

**Example:**

.. code-block:: python

    engine.set_species_lower_bound(8, 400, "moles")
)doc")
            .def("set_species_upper_bound", py::overload_cast<const std::string&, double, const std::string&, std::optional<std::string>>
                 (&ChemicalEngineMaps::set_species_upper_bound),
             py::arg("species"), py::arg("val"), py::arg("units")="moles", py::arg("phase") = py::none(),
             R"doc(
Sets an upper bound for a species identified by name.

:param str species: Species name.
:param float val: Upper limit in units.
:param str units: Units of amount ("moles", "kg", "m3"), default "moles".
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

    engine.set_species_upper_bound('CaOH+', 500, "kg")
)doc")
            .def("set_species_upper_bound", py::overload_cast<Index, double, const std::string&>
                 (&ChemicalEngineMaps::set_species_upper_bound),
             py::arg("ispecies"), py::arg("val"), py::arg("units")="moles",
             R"doc(
Sets an upper bound for a species identified by its index.

:param int ispecies: Index of the species.
:param float val: Upper limit in units.
:param str units: Units of amount ("moles", "kg", "m3"), default "moles".

**Example:**

.. code-block:: python

    engine.set_species_upper_bound(8, 400, "moles")
)doc")

            .def("suppress_phase", &ChemicalEngineMaps::suppress_phase,
             py::arg("phase_name"), py::arg("min_amount")=0, py::arg("max_amount")=1e-15,
             R"doc(
Suppresses a phase in GEM calculation.

:param str phase_name: Phase name.
:param float min_amount: Lower amount of species in mol, default 0.
:param float max_amount: Upper amount of species in mol, default 1e-15.

**Example:**

.. code-block:: python

    engine.suppress_phase('gas_gen')
)doc")
            .def("suppress_multiple_phases", &ChemicalEngineMaps::suppress_multiple_phases,
             py::arg("phase_name_list"), py::arg("min_amount")=0, py::arg("max_amount")=1e-15,
             R"doc(
Suppresses multiple phases in calculation as given in phase names list.

:param list phase_name_list: Phase name list.
:param float min_amount: Lower amount of species in mol, default 0.
:param float max_amount: Upper amount of species in mol, default 1e-15.

**Example:**

.. code-block:: python

    engine.suppress_multiple_phases(['Dolomite-dis', 'Tin'])
)doc")
            .def("suppress_species", &ChemicalEngineMaps::suppress_species, py::arg("species_name"),
             py::arg("min_amount")=0, py::arg("max_amount")=1e-15, py::arg("phase") = py::none(),
             R"doc(
Suppresses a species in calculation.

:param str species_name: Species name.
:param float min_amount: Lower amount of species in mol, default 0.
:param float max_amount: Upper amount of species in mol, default 1e-15.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

    engine.suppress_species('Ca(CO3)@')
)doc")
            .def("suppress_multiple_species", &ChemicalEngineMaps::suppress_multiple_species,
             py::arg("species_list"), py::arg("min_amount")=0, py::arg("max_amount")=1e-15,
             R"doc(
Suppresses multiple species in GEM calculation as given in species name list.

:param list species_list: Species name list.
:param float min_amount: Lower amount of species in mol, default 0.
:param float max_amount: Upper amount of species in mol, default 1e-15.

**Example:**

.. code-block:: python

    engine.suppress_multiple_species(['ClO4-', 'Cl-'])
)doc")
            .def("activate_phase", &ChemicalEngineMaps::activate_phase,
             R"doc(
Activate a suppressed phase in GEM calculation.

:param str phase_name: Phase name.

**Example:**

.. code-block:: python

    engine.activate_phase('gas_gen')
)doc")
            .def("activate_multiple_phases", &ChemicalEngineMaps::activate_multiple_phases,
             R"doc(
Activate multiple suppressed phases given in list.

:param list phase_name_list: Phase name list.

**Example:**

.. code-block:: python

    engine.activate_multiple_phases(['Dolomite-dis', 'Tin'])
)doc")
            .def("activate_species", &ChemicalEngineMaps::activate_species,
             py::arg("species_name"), py::arg("phase") = py::none(),
             R"doc(
Activate a suppressed species in phase.

:param str species_name: Species name.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

    engine.activate_species('Ca(CO3)@')
)doc")
            .def("activate_multiple_species", &ChemicalEngineMaps::activate_multiple_species,
             R"doc(
Activate multiple suppressed species given in the list.

:param list species_list: Species name list.

**Example:**

.. code-block:: python

    engine.activate_multiple_species(['ClO4-', 'Cl-'])
)doc")
            .def("__repr__", [](const ChemicalEngineMaps &self)
             {
                std::stringstream ss;
                ss << self;
                return ss.str();
             }, R"doc(
Returns a string representation of the ChemicalEngineDicts instance.

Reflects the current staged T, P, and bulk composition even before
``equilibrate()`` is called.

**Example:**

.. code-block:: python

    print(engine)
)doc")
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
