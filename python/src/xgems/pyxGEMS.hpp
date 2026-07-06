// xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2018-2025 Allan Leal, Dmtrii Kulik, G.D. Miron, S.Dmytriieva
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

void exportChemicalEngine(py::module &m)
{
    auto reequilibrate1 = static_cast<int (ChemicalEngine::*)()>(&ChemicalEngine::reequilibrate);
    auto reequilibrate2 = static_cast<int (ChemicalEngine::*)(bool)>(&ChemicalEngine::reequilibrate);

    // Bind ChemicalEngineOptions
    py::class_<ChemicalEngineOptions>(m, "ChemicalEngineOptions",
        R"doc(
Configuration options for the ChemicalEngine.

Attributes:
    warmstart (bool): Use smart start for faster convergence.
    print_zero_amounts (bool): Whether to print zero-amount species/phases.
)doc")
        .def(py::init<>())
        .def_readwrite("warmstart", &ChemicalEngineOptions::warmstart,
            R"doc(
Whether to use smart start for faster convergence (may be less accurate).

:type: bool
)doc")
        .def_readwrite("print_zero_amounts", &ChemicalEngineOptions::print_zero_amounts,
            R"doc(
Whether to include zero-amount species/phases in the printed output.

:type: bool
)doc");

    m.def("update_loggers", &update_loggers,
          R"doc(
  Updates logger settings for the ChemicalEngine.
  
  :param bool use_cout: If true, prints log messages to stdout.
  :param str logfile_name: The filename for a rotating log file. If empty, file logging is disabled.
  :param int log_level: Verbosity level for logging (trace=0, debug=1, info=2, warn=3, err=4, critical=5, off=6).
  
  **Example:**
  
  .. code-block:: python
  
      update_loggers(True, "chemical_log.txt", 2)
  )doc",
          py::arg("use_cout"), py::arg("logfile_name"), py::arg("log_level"));

    py::class_<ChemicalEngine>(m, "ChemicalEngine")
        .def(py::init<>(),
             R"doc(
  Default constructor for the ChemicalEngine.
  
  Creates an empty ChemicalEngine instance. Use one of the initialization routines
  (such as ``initialize()`` or ``initializeFromJsonStrings()``) to configure the system.
  
  **Example:**
  
  .. code-block:: python
  
      engine = ChemicalEngine()
      engine.initialize("my-system-dat.lst")
  )doc")
        .def(py::init<std::string>(), py::arg("filename"),
             R"doc(
  Constructs a ChemicalEngine instance by loading a GEM-Selektor project file.
  
  :param str filename: The file path for the chemical system definition (e.g., "my-system-dat.lst").
  
  **Example:**
  
  .. code-block:: python
  
      engine = ChemicalEngine("my-system-dat.lst")
  )doc")
        .def("initialize", &ChemicalEngine::initialize, py::arg("filename"),
             R"doc(
  Initializes the ChemicalEngine from a GEM-Selektor project file.
  
  :param str filename: Path to the system definition file (e.g., "my-system-dat.lst").
  
  **Example:**
  
  .. code-block:: python
  
      engine.initialize("my-system-dat.lst")
  )doc")
        .def("initializeFromJsonStrings", &ChemicalEngine::initializeFromJsonStrings,
             py::arg("dch_json"), py::arg("ipm_json"), py::arg("dbr_json"),
             R"doc(
  Initializes the ChemicalEngine using JSON strings.
  
  :param str dch_json: JSON string describing the chemical system.
  :param str ipm_json: JSON string with parameter and algorithm settings.
  :param str dbr_json: JSON string containing node composition details.
  
  **Example:**
  
  .. code-block:: python
  
      engine.initializeFromJsonStrings(dch, ipm, dbr)
  )doc")
        .def("readDbrFromFile", &ChemicalEngine::readDbrFromFile, py::arg("filename"),
             R"doc(
  Reads the input node composition (DBR) from a file.
  
  :param str filename: Path to the DBR file.
  
  **Example:**
  
  .. code-block:: python
  
      engine.readDbrFromFile("input.dbr")
  )doc")
        .def("readDbrFromJsonString", &ChemicalEngine::readDbrFromJsonString, py::arg("dbr_json"),
             R"doc(
  Reads the input node composition (DBR) from a JSON string.
  
  :param str dbr_json: JSON string containing DBR composition data.
  
  **Example:**
  
  .. code-block:: python
  
      engine.readDbrFromJsonString("{}")
  )doc")
        .def("writeDbrToFile", &ChemicalEngine::writeDbrToFile, py::arg("filename"),
             R"doc(
  Writes the current node composition (DBR) to a file.
  
  :param str filename: Path to the output DBR file.
  
  **Example:**
  
  .. code-block:: python
  
      engine.writeDbrToFile("output.dbr")
  )doc")
        .def("writeDbrToJsonString", &ChemicalEngine::writeDbrToJsonString,
             R"doc(
  Returns the current DBR as a JSON string.
  
  :return: JSON string representing the system composition.
  :rtype: str
  
  **Example:**
  
  .. code-block:: python
  
      json_state = engine.writeDbrToJsonString()
      print(json_state)
  )doc")
        .def("numElements", &ChemicalEngine::numElements,
             R"doc(
  Returns the number of elements in the chemical system.
  
  **Example:**
  
  .. code-block:: python
  
      num_elements = engine.numElements()
      print(num_elements)
  )doc")
        .def("numSpecies", &ChemicalEngine::numSpecies,
             R"doc(
  Returns the number of species in the chemical system.
  
  **Example:**
  
  .. code-block:: python
  
      num_species = engine.numSpecies()
      print(num_species)
  )doc")
        .def("numPhases", &ChemicalEngine::numPhases,
             R"doc(
  Returns the number of phases in the chemical system.
  
  **Example:**
  
  .. code-block:: python
  
      num_phases = engine.numPhases()
      print(num_phases)
  )doc")
        .def("numSpeciesInPhase", static_cast<Index (ChemicalEngine::*)(Index) const>(&ChemicalEngine::numSpeciesInPhase),
             py::arg("iphase"),
             R"doc(
  Returns the number of species in a given phase.
  Access specified element with bounds checking.
  
  :param int iphase: Index of the phase.
  
  **Example:**
  
  .. code-block:: python
  
      num_species_in_phase = engine.numSpeciesInPhase(0)
      print(num_species_in_phase)
  )doc")
        .def("numSpeciesInPhase", static_cast<Index (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::numSpeciesInPhase),
             py::arg("phase"),
             R"doc(
  Returns the number of species in in a phase identified by name.

  :param str phase: Name of the phase.

  **Example:**

  .. code-block:: python

      num_species_in_phase = engine.numSpeciesInPhase("aq_gen")
      print(num_species_in_phase)
  )doc")
        .def("numSpeciesInPhase_i", &ChemicalEngine::numSpeciesInPhase_i,
             py::arg("iphase"),
             R"doc(
  Returns the number of species in a given phase.
  Access specified element without bounds checking.

  :param int iphase: Index of the phase.

  **Example:**

  .. code-block:: python

      num_species_in_phase = engine.numSpeciesInPhase_i(0)
      print(num_species_in_phase)
  )doc")
        .def("elementName", &ChemicalEngine::elementName, py::arg("index"),
             R"doc(
  Returns the name of an element by its index.
  Access specified element with bounds checking.

  :param int index: Index of the element.

  **Example:**

  .. code-block:: python

      element_name = engine.elementName(0)
      print(element_name)
  )doc")
        .def("elementName_i", &ChemicalEngine::elementName_i, py::arg("index"),
             R"doc(
  Returns the name of an element by its index.
  Access specified element without bounds checking.

  :param int index: Index of the element.

  **Example:**

  .. code-block:: python

      element_name = engine.elementName_i(0)
      print(element_name)
  )doc")
        .def("speciesName", &ChemicalEngine::speciesName, py::arg("index"),
             R"doc(
  Returns the name of a species by its index.
  Access specified element with bounds checking.
  
  :param int index: Index of the species.
  
  **Example:**
  
  .. code-block:: python
  
      species_name = engine.speciesName(0)
      print(species_name)
  )doc")
        .def("speciesName_i", &ChemicalEngine::speciesName_i, py::arg("index"),
             R"doc(
  Returns the name of a species by its index.
  Access specified element without bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      species_name = engine.speciesName_i(0)
      print(species_name)
  )doc")
        .def("speciesCharge", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::speciesCharge), py::arg("index"),
             R"doc(
  Returns the electrical charge of a species by its index.
  Access specified element with bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      charge = engine.speciesCharge(0)
      print(charge)
  )doc")
        .def("speciesCharge", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::speciesCharge),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
  Returns the electrical charge of a species identified by name.

  :param str species: Name of the species.
  :param str phase: Name of the phase the species was included in. If None get the first index.

  **Example:**

  .. code-block:: python

      charge = engine.speciesCharge("SiO2")
      print(charge)
  )doc")
        .def("speciesCharge_i", &ChemicalEngine::speciesCharge_i, py::arg("index"),
             R"doc(
  Returns the electrical charge of a species by its index.
  Access specified element without bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      charge = engine.speciesCharge_i(0)
      print(charge)
  )doc")
        .def("phaseName", &ChemicalEngine::phaseName, py::arg("index"),
             R"doc(
  Returns the name of a phase by its index.
  Access specified element with bounds checking.

  :param int index: Index of the phase.

  **Example:**

  .. code-block:: python

      phase_name = engine.phaseName(0)
      print(phase_name)
  )doc")
        .def("phaseName_i", &ChemicalEngine::phaseName_i, py::arg("index"),
             R"doc(
  Returns the name of a phase by its index.
  Access specified element without bounds checking.

  :param int index: Index of the phase.

  **Example:**

  .. code-block:: python

      phase_name = engine.phaseName_i(0)
      print(phase_name)
  )doc")
        .def("indexElement", &ChemicalEngine::indexElement, py::arg("name"),
             R"doc(
  Returns the index of an element by its name.
  
  :param str name: Name of the element.
  
  **Example:**
  
  .. code-block:: python
  
      index = engine.indexElement("O")
      print(index)
  )doc")
        .def("indexSpecies", &ChemicalEngine::indexSpecies, py::arg("name"), py::arg("phase") = py::none(),
             R"doc(
Returns the index of a species by its name.

:param str name: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  index = engine.indexSpecies("H2O@")
  print(index)
)doc")

        .def("indexSpeciesAll", &ChemicalEngine::indexSpeciesAll, py::arg("name"),
             R"doc(
Returns all indices of species matching the specified name.

:param str name: Name of the species.

**Example:**

.. code-block:: python

  indices = engine.indexSpeciesAll("H2O@")
  print(indices)
)doc")
        .def("indexSpeciesMap", &ChemicalEngine::indexSpeciesMap, py::arg("name"),
             R"doc(
Returns dictionary of all phase name and corresponding species index.

:param str name: Name of the species.

**Example:**

.. code-block:: python

  dict = engine.indexSpeciesMap("Anorthite")
  print(dict)
)doc")
        .def("indexPhase", &ChemicalEngine::indexPhase, py::arg("name"),
             R"doc(
Returns the index of a phase by its name.

:param str name: Name of the phase.

**Example:**

.. code-block:: python

  index = engine.indexPhase("aq_gen")
  print(index)
)doc")

        .def("indexPhaseAll", &ChemicalEngine::indexPhaseAll, py::arg("name"),
             R"doc(
Returns all indices of phases matching the specified name.

:param str name: Name of the phase.

**Example:**

.. code-block:: python

  indices = engine.indexPhaseAll("aq_gen")
  print(indices)
)doc")

        .def("indexPhaseWithSpecies", &ChemicalEngine::indexPhaseWithSpecies, py::arg("species_index"),
             R"doc(
Returns the index of the phase containing a given species.

:param int species_index: Index of the species.

**Example:**

.. code-block:: python

  phase_index = engine.indexPhaseWithSpecies(0)
  print(phase_index)
)doc")

        .def("indexFirstSpeciesInPhase", static_cast<Index (ChemicalEngine::*)(Index) const>(&ChemicalEngine::indexFirstSpeciesInPhase),
             py::arg("iphase"),
             R"doc(
Returns the index of the first species in a specified phase.
Access specified element with bounds checking.

:param int iphase: Index of the phase.

**Example:**

.. code-block:: python

  species_index = engine.indexFirstSpeciesInPhase(0)
  print(species_index)
)doc")
        .def("indexFirstSpeciesInPhase", static_cast<Index (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::indexFirstSpeciesInPhase),
             py::arg("phase"),
             R"doc(
Returns the index of the first species in a specified phase.

:param str phase: Name of the phase.

**Example:**

.. code-block:: python

  species_index = engine.indexFirstSpeciesInPhase("aq_gen")
  print(species_index)
)doc")
        .def("indexFirstSpeciesInPhase_i", &ChemicalEngine::indexFirstSpeciesInPhase_i, py::arg("iphase"),
             R"doc(
Returns the index of the first species in a specified phase.
Access specified element without bounds checking.

:param int iphase: Index of the phase.

**Example:**

.. code-block:: python

  species_index = engine.indexFirstSpeciesInPhase_i(0)
  print(species_index)
)doc")

        .def("elementMolarMasses", &ChemicalEngine::elementMolarMasses,
             R"doc(
Returns the molar masses of all elements in the system (kg/mol).

**Example:**

.. code-block:: python

  molar_masses = engine.elementMolarMasses()
  print(molar_masses)
)doc",
             py::return_value_policy::reference_internal)

        .def("speciesMolarMasses", &ChemicalEngine::speciesMolarMasses,
             R"doc(
Returns the molar masses of all species in the system (kg/mol).

**Example:**

.. code-block:: python

  species_molar_masses = engine.speciesMolarMasses()
  print(species_molar_masses)
)doc",
             py::return_value_policy::reference_internal)

        .def("formulaMatrix", &ChemicalEngine::formulaMatrix,
             R"doc(
Returns the formula matrix of the system.

The formula matrix represents the stoichiometric coefficients of elements in species.

**Example:**

.. code-block:: python

  formula_matrix = engine.formulaMatrix()
  print(formula_matrix)
)doc",
             py::return_value_policy::reference_internal)
             .def_property("options",
                &ChemicalEngine::options,
                &ChemicalEngine::setOptions,
                R"doc(
            Gets or sets the options used by the ChemicalEngine.
            
            :type: ChemicalEngineOptions
            
            **Example (get):**
            
            .. code-block:: python
            
              opts = engine.options
              print("Warmstart:", opts.warmstart)
            
            **Example (set):**
            
            .. code-block:: python
            
              opts = engine.options
              opts.warmstart = False
              opts.print_zero_amounts = True # prints phases and species with zero amounts
              engine.options = opts
            )doc")

        .def("setWarmStart", &ChemicalEngine::setWarmStart,
             R"doc(
Enables warm start for the ChemicalEngine.

Warm start uses the previous solution as an initial guess for faster convergence.

**Example:**

.. code-block:: python

  engine.setWarmStart()
)doc")

        .def("setColdStart", &ChemicalEngine::setColdStart,
             R"doc(
Enables cold start for the ChemicalEngine.

Cold start resets the initial guess to default values for robust convergence.

**Example:**

.. code-block:: python

  engine.setColdStart()
)doc")

        .def("setSpeciesUpperLimit", static_cast<void (ChemicalEngine::*)(std::string, double, std::optional<std::string>)>(&ChemicalEngine::setSpeciesUpperLimit),
             py::arg("name"), py::arg("limit"), py::arg("phase") = py::none(),
             R"doc(
Sets the upper limit for a species by its name.

:param str name: Name of the species.
:param float limit: Upper limit for the species amount.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setSpeciesUpperLimit("H2O@", 10.0)
)doc")
        .def("setSpeciesUpperLimit", static_cast<void (ChemicalEngine::*)(Index, double)>(&ChemicalEngine::setSpeciesUpperLimit),
             py::arg("index"), py::arg("limit"),
             R"doc(
Sets the upper limit for a species by its index.
Access specified element with bounds checking.

:param int index: Index of the species.
:param float limit: Upper limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesUpperLimit(0, 10.0)
)doc")
        .def("setSpeciesUpperLimit_i", &ChemicalEngine::setSpeciesUpperLimit_i, py::arg("index"), py::arg("limit"),
             R"doc(
Sets the upper limit for a species by its index.
Access specified element without bounds checking.


:param int index: Index of the species.
:param float limit: Upper limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesUpperLimit_i(0, 10.0)
)doc")
        .def("setSpeciesLowerLimit", static_cast<void (ChemicalEngine::*)(std::string, double, std::optional<std::string>)>(&ChemicalEngine::setSpeciesLowerLimit),
             py::arg("name"), py::arg("limit"), py::arg("phase") = py::none(),
             R"doc(
Sets the lower limit for a species by its name.

:param str name: Name of the species.
:param float limit: Lower limit for the species amount.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setSpeciesLowerLimit("SiO2", 0.1)
)doc")
        .def("setSpeciesLowerLimit", static_cast<void (ChemicalEngine::*)(Index, double)>(&ChemicalEngine::setSpeciesLowerLimit),
             py::arg("index"), py::arg("limit"),
             R"doc(
Sets the lower limit for a species by its index.
Access specified element with bounds checking.

:param int index: Index of the species.
:param float limit: Lower limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesLowerLimit(0, 0.1)
)doc")
        .def("setSpeciesLowerLimit_i", &ChemicalEngine::setSpeciesLowerLimit_i,
             py::arg("index"), py::arg("limit"),
             R"doc(
Sets the lower limit for a species by its index.
Access specified element without bounds checking.

:param int index: Index of the species.
:param float limit: Lower limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesLowerLimit_i(0, 0.1)
)doc")
        .def("setSpeciesAmount", static_cast<void (ChemicalEngine::*)(std::string, double, std::optional<std::string>)>(&ChemicalEngine::setSpeciesAmount),
             py::arg("name"), py::arg("amount"), py::arg("phase") = py::none(),
             R"doc(
Sets the amount of a species by its name.

:param str name: Name of the species.
:param float amount: Amount of the species in moles.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setSpeciesAmount("SiO2", 1.0)
)doc")
        .def("setSpeciesAmount", static_cast<void (ChemicalEngine::*)(Index, double)>(&ChemicalEngine::setSpeciesAmount),
             py::arg("ispecies"), py::arg("amount"),
             R"doc(
Sets the amount of a species by its index.
Access specified element with bounds checking.

:param int ispecies: Index of the species.
:param float amount: Amount of the species in moles.

**Example:**

.. code-block:: python

  engine.setSpeciesAmount(0, 1.0)
)doc")
        .def("setSpeciesAmount_i", &ChemicalEngine::setSpeciesAmount_i, py::arg("ispecies"), py::arg("amount"),
             R"doc(
Sets the amount of a species by its index.
Access specified element without bounds checking.

:param int ispecies: Index of the species.
:param float amount: Amount of the species in moles.

**Example:**

.. code-block:: python

  engine.setSpeciesAmount_i(0, 1.0)
)doc")

        .def("setStandardMolarGibbsEnergy", static_cast<void (ChemicalEngine::*)(std::string, double, std::optional<std::string>)>(&ChemicalEngine::setStandardMolarGibbsEnergy),
py::arg("name"), py::arg("value"), py::arg("phase") = py::none(),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).

:param str name: Name of the species.
:param float value: Standard molar Gibbs energy value (J/mol).
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setStandardMolarGibbsEnergy("H2O", -237.13)
)doc")
        .def("setStandardMolarGibbsEnergy", static_cast<void (ChemicalEngine::*)(Index, double)>(&ChemicalEngine::setStandardMolarGibbsEnergy),
             py::arg("ispecies"), py::arg("value"),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).
Access specified element with bounds checking.

:param int ispecies: Index of the species.
:param float value: Standard molar Gibbs energy value (J/mol).

**Example:**

.. code-block:: python

  engine.setStandardMolarGibbsEnergy(0, -237.13)
)doc")
        .def("setStandardMolarGibbsEnergy_i", &ChemicalEngine::setStandardMolarGibbsEnergy_i,
             py::arg("ispecies"), py::arg("value"),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).
Access specified element with bounds checking.

:param int ispecies: Index of the species.
:param float value: Standard molar Gibbs energy value (J/mol).

**Example:**

.. code-block:: python

  engine.setStandardMolarGibbsEnergy_i(0, -237.13)
)doc")
    .def("setPT", &ChemicalEngine::setPT,
             R"doc(
                   Sets the temperature and pressure of the system.
                   
                   :param float temperature: Temperature in Kelvin (K).
                   :param float pressure: Pressure in Pascals (Pa).

                    :return bool: True if PT was set correctly, False if out of range.
                   
                   **Example:**
                   
                   .. code-block:: python
                   
                       engine.setPT(298.15, 101325)
              )doc",
             py::arg("temperature"), py::arg("pressure"))

    .def("setB", &ChemicalEngine::setB,
             R"doc(
                   Sets the bulk composition of the system.
                   
                   :param list[float] composition: List of element amounts (mol).
                   
                   Set element amounts of a solution with 0.001 M MgCl2 and 0.001 M CO2:
                   - **elements order:** C, Ca, Cl, H, Mg, O, Zz (charge)
                   - **Note:** The order of elements in the list should match the order in the chemical system.
                   
                   **Example:**
                   
                   .. code-block:: python
                   
                       engine.setB([0.001, 1e-09, 0.004, 110.6837, 0.002, 55.34405, 0])
              )doc",
             py::arg("composition"))

    .def("reequilibrate", reequilibrate1,
             R"doc(
                   Re-equilibrates the system using the current state as the initial guess. Default with automatic inital guess.
                   
                   **Example:**
                   
                   .. code-block:: python
                   
                       engine.reequilibrate()
                  
                   **Return Codes**

                   The function returns an integer code indicating the status:
                   
                   - 0: No GEM re-calculation needed
                   - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
                   - 2: OK after GEM calculation with LPP AIA
                   - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
                   - 4: Failure (no result) in GEM calculation with LPP AIA
                   - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
                   - 6: OK after GEM calculation with SIA
                   - 7: Bad (not fully trustful) result after GEM calculation with SIA
                   - 8: Failure (no result) in GEM calculation with SIA
                   - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.

              )doc")

    .def("reequilibrate", reequilibrate2,
             R"doc(
                   Re-equilibrates the system using the current state as the initial guess.
   
                   :param bool warmstart: equilibration with smart initial guess.
                   
                   **Example:**
                   
                   .. code-block:: python
                   
                       engine.reequilibrate(True)
                  
                   **Return Codes**

                   The function returns an integer code indicating the status:
                   
                   - 0: No GEM re-calculation needed
                   - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
                   - 2: OK after GEM calculation with LPP AIA
                   - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
                   - 4: Failure (no result) in GEM calculation with LPP AIA
                   - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
                   - 6: OK after GEM calculation with SIA
                   - 7: Bad (not fully trustful) result after GEM calculation with SIA
                   - 8: Failure (no result) in GEM calculation with SIA
                   - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.

              )doc")

    .def("equilibrate", &ChemicalEngine::equilibrate,
                R"doc(
           Compute the chemical equilibrium state from temperature, pressure, and element amounts.
           
           This method resets any existing state and computes a new equilibrium based on the
           provided temperature (in Kelvin), pressure (in Pascals), and vector of element amounts (in mol).
           The order of the elements in the vector must match the chemical system's element ordering.
           
           :param float T: Temperature in Kelvin.
           :param float P: Pressure in Pascals.
           :param b: Vector of element amounts (in mol). Must be a NumPy array or Eigen vector of correct length.
           :type b: numpy.ndarray or Eigen::VectorXd
           :returns: Integer return code indicating the status of the equilibrium computation.
           :rtype: int
           
           **Example**
           
           .. code-block:: python
           
               import numpy as np
           
               # Element amounts in mol (order: C, Ca, Cl, H, Mg, O, Zz)
               b = np.array([0.001, 1e-9, 0.004, 110.6837, 0.002, 55.34405, 0.0])
           
               # Compute equilibrium
               ret = engine.equilibrate(298.15, 101325, b)
           
           **Return Codes**
           
           The function returns an integer code indicating the status:
           
           - 0: No GEM re-calculation needed
           - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
           - 2: OK after GEM calculation with LPP AIA
           - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
           - 4: Failure (no result) in GEM calculation with LPP AIA
           - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
           - 6: OK after GEM calculation with SIA
           - 7: Bad (not fully trustful) result after GEM calculation with SIA
           - 8: Failure (no result) in GEM calculation with SIA
           - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.
           
           .. note::
              Ensure the `b` vector is properly ordered and matches the number of elements
              defined in the chemical system.
           )doc")
           
           
        .def("converged", &ChemicalEngine::converged,
             R"doc(
            Checks if the equilibrium computation has converged.
            
            :return bool: True if the computation has converged, False otherwise.
            
            **Example:**
            
            .. code-block:: python
            
                is_converged = engine.converged()
                print(is_converged)
            )doc")

        .def("numIterations", &ChemicalEngine::numIterations,
             R"doc(
            Returns the number of iterations performed during the last equilibrium computation.
            
            :return int: Number of iterations.
            
            **Example:**
            
            .. code-block:: python
            
                iterations = engine.numIterations()
                print(iterations)
            )doc")

        .def("elapsedTime", &ChemicalEngine::elapsedTime,
             R"doc(
            Returns the elapsed time for the last equilibrium computation (in seconds).
            
            :return float: Elapsed time in seconds.
            
            **Example:**
            
            .. code-block:: python
            
                time = engine.elapsedTime()
                print(time)
            )doc")

        .def("temperature", &ChemicalEngine::temperature,
             R"doc(
            Returns the current temperature of the system (K).
            
            :return float: Temperature in Kelvin.
            
            **Example:**
            
            .. code-block:: python
            
                temp = engine.temperature()
                print(temp)
            )doc")

        .def("pressure", &ChemicalEngine::pressure,
             R"doc(
            Returns the current pressure of the system (Pa).
            
            :return float: Pressure in Pascals.
            
            **Example:**
            
            .. code-block:: python
            
                pressure = engine.pressure()
                print(pressure)
            )doc")

        .def("elementAmounts", &ChemicalEngine::elementAmounts,
             R"doc(
            Returns the amounts of all elements in the system (mol).
            
            **Example:**
            
            .. code-block:: python
            
                amounts = engine.elementAmounts()
                print(amounts)
            )doc",
             py::return_value_policy::reference_internal)

        .def("elementAmountsInPhase", static_cast<Vector (ChemicalEngine::*)(Index) const>(&ChemicalEngine::elementAmountsInPhase),
             py::arg("iphase"),
             R"doc(
            Returns the amounts of all elements in a specific phase (mol).
            Access specified element with bounds checking.
            
            :param int phase_index: Index of the phase.
            
            **Example:**
            
            .. code-block:: python
            
                amounts = engine.elementAmountsInPhase(0)
                print(amounts)
            )doc")
        .def("elementAmountsInPhase", static_cast<Vector (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::elementAmountsInPhase),
             py::arg("phase"),
             R"doc(
            Returns the amounts of all elements in a specific phase (mol).

            :param str phase: Name of the phase.

            **Example:**

            .. code-block:: python

                amounts = engine.elementAmountsInPhase("aq_gen"")
                print(amounts)
            )doc")
        .def("elementAmountsInPhase_i", &ChemicalEngine::elementAmountsInPhase_i,
             py::arg("iphase"),
             R"doc(
            Returns the amounts of all elements in a specific phase (mol).
            Access specified element without bounds checking.

            :param int iphase: Index of the phase.

            **Example:**

            .. code-block:: python

                amounts = engine.elementAmountsInPhase_i(0)
                print(amounts)
            )doc")

        .def("elementAmountsInSpecies", &ChemicalEngine::elementAmountsInSpecies,
             R"doc(
            Returns the amounts of all elements in a specific species (mol).
            
            :param int species_index: Index of the species.
            
            **Example:**
            
            .. code-block:: python
            
                amounts = engine.elementAmountsInSpecies(0)
                print(amounts)
            )doc",
             py::arg("species_index"))

        .def("speciesAmount", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::speciesAmount),
             py::arg("ispecies"),
             R"doc(
            Returns the amount of a species by its index (mol).
            Access specified element with bounds checking.
            
            :param int ispecies: Index of the species.
            
            **Example:**
            
            .. code-block:: python
            
                amount = engine.speciesAmount(0)
                print(amount)
            )doc")
        .def("speciesAmount", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::speciesAmount),
            py::arg("name"), py::arg("phase") = py::none(),
             R"doc(
            Returns the amount of a species by its name (mol).
            
            :param str name: Name of the species.
            :param str phase: Name of the phase the species was included in. If None get the first index.
            
            **Example:**
            
            .. code-block:: python
            
                amount = engine.speciesAmount("H2O@")
                print(amount)
            )doc")
        .def("speciesAmount_i", &ChemicalEngine::speciesAmount_i, py::arg("ispecies"),
             R"doc(
            Returns the amount of a species by its index (mol).
            Access specified element without bounds checking.

            :param int ispecies: Index of the species.

            **Example:**

            .. code-block:: python

                amount = engine.speciesAmount_i(0)
                print(amount)
            )doc")

        .def("setSpeciesAmounts", &ChemicalEngine::setSpeciesAmounts,
             R"doc(
            Sets the amounts of all species in the system (mol).
            
            :param list[float] amounts: List of species amounts.
            
            **Example:**
            
            .. code-block:: python
            
                engine.setSpeciesAmounts([1.0, 0.5, 0.2])
            )doc",
             py::arg("amounts"), py::return_value_policy::reference_internal)

        .def("speciesAmounts", &ChemicalEngine::speciesAmounts,
             R"doc(
            Returns the amounts of all species in the system (mol).
            
            **Example:**
            
            .. code-block:: python
            
                amounts = engine.speciesAmounts()
                print(amounts)
            )doc",
             py::return_value_policy::reference_internal)

        .def("setSpeciesUpperLimits", &ChemicalEngine::setSpeciesUpperLimits,
             R"doc(
            Sets the upper limits for all species in the system (mol).
            
            :param list[float] limits: List of upper limits for species amounts.
            
            **Example:**
            
            .. code-block:: python
            
                engine.setSpeciesUpperLimits([10.0, 5.0, 2.0])
            )doc",
             py::arg("limits"), py::return_value_policy::reference_internal)

        .def("setSpeciesLowerLimits", &ChemicalEngine::setSpeciesLowerLimits,
             R"doc(
            Sets the lower limits for all species in the system (mol).
            
            :param list[float] limits: List of lower limits for species amounts.
            
            **Example:**
            
            .. code-block:: python
            
                engine.setSpeciesLowerLimits([0.1, 0.05, 0.02])
            )doc",
             py::arg("limits"), py::return_value_policy::reference_internal)

        .def("speciesUpperLimits", &ChemicalEngine::speciesUpperLimits,
             R"doc(
            Returns the upper limits for all species in the system (mol).
            
            **Example:**
            
            .. code-block:: python
            
                limits = engine.speciesUpperLimits()
                print(limits)
            )doc",
             py::return_value_policy::reference_internal)

        .def("speciesLowerLimits", &ChemicalEngine::speciesLowerLimits,
             R"doc(
            Returns the lower limits for all species in the system (mol).
            
            **Example:**
            
            .. code-block:: python
            
                limits = engine.speciesLowerLimits()
                print(limits)
            )doc",
             py::return_value_policy::reference_internal)

        .def("moleFractions", &ChemicalEngine::moleFractions,
             R"doc(
            Returns the mole fractions of all species in the system.
            
            **Example:**
            
            .. code-block:: python
            
                fractions = engine.moleFractions()
                print(fractions)
            )doc",
             py::return_value_policy::reference_internal)

        .def("speciesMolalities", &ChemicalEngine::speciesMolalities,
             R"doc(
            Returns the molalities of all species in the system (mol/kg).
            
            **Example:**
            
            .. code-block:: python
            
                molalities = engine.speciesMolalities()
                print(molalities)
            )doc",
             py::return_value_policy::reference_internal)

        .def("lnActivityCoefficients", &ChemicalEngine::lnActivityCoefficients,
             R"doc(
            Returns the natural logarithms of activity coefficients for all species.
            
            **Example:**
            
            .. code-block:: python
            
                ln_coeffs = engine.lnActivityCoefficients()
                print(ln_coeffs)
            )doc",
             py::return_value_policy::reference_internal)

        .def("lnActivities", &ChemicalEngine::lnActivities,
             R"doc(
            Returns the natural logarithms of activities for all species.
            
            **Example:**
            
            .. code-block:: python
            
                ln_activities = engine.lnActivities()
                print(ln_activities)
            )doc",
             py::return_value_policy::reference_internal)

        .def("lnConcentrations", &ChemicalEngine::lnConcentrations,
             R"doc(
            Returns the natural logarithms of concentrations for all species.
            
            **Example:**
            
            .. code-block:: python
            
                ln_concentrations = engine.lnConcentrations()
                print(ln_concentrations)
            )doc",
             py::return_value_policy::reference_internal)

        .def("chemicalPotentials", &ChemicalEngine::chemicalPotentials,
             R"doc(
                Returns the chemical potentials of all species in the system (J/mol).
                
                **Example:**
                
                .. code-block:: python
                
                    potentials = engine.chemicalPotentials()
                    print(potentials)
                )doc",
             py::return_value_policy::reference_internal)

        .def("standardMolarGibbsEnergy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::standardMolarGibbsEnergy),
            py::arg("ispecies"),
             R"doc(
                Returns the standard molar Gibbs energy of a specific species by its index (J/mol).
                Access specified element with bounds checking.
                
                :param int ispecies: Index of the species.
                
                **Example:**
                
                .. code-block:: python
                
                    gibbs_energy = engine.standardMolarGibbsEnergy(0)
                    print(gibbs_energy)
                )doc")
        .def("standardMolarGibbsEnergy", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::standardMolarGibbsEnergy),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
                Returns the standard molar Gibbs energy of a specific species by its index (J/mol).

                :param str species: Name of the species.
                :param str phase: Name of the phase the species was included in. If None get the first index.

                **Example:**

                .. code-block:: python

                    gibbs_energy = engine.standardMolarGibbsEnergy("SiO2")
                    print(gibbs_energy)
                )doc")
        .def("standardMolarGibbsEnergy_i", &ChemicalEngine::standardMolarGibbsEnergy_i,
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar Gibbs energy of a specific species by its index (J/mol).
                Access specified element without bounds checking.

                :param int ispecies: Index of the species.

                **Example:**

                .. code-block:: python

                    gibbs_energy = engine.standardMolarGibbsEnergy_i(0)
                    print(gibbs_energy)
                )doc")

        .def("standardMolarEnthalpy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::standardMolarEnthalpy),
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar enthalpy of a specific species by its index (J/mol).
                Access specified element with bounds checking.
                
                :param int ispecies: Index of the species.
                
                **Example:**
                
                .. code-block:: python
                
                    enthalpy = engine.standardMolarEnthalpy(0)
                    print(enthalpy)
                )doc")
        .def("standardMolarEnthalpy", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::standardMolarEnthalpy),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
                Returns the standard molar enthalpy of a specific species by its index (J/mol).

                :param str species: Name of the species.
                :param str phase: Name of the phase the species was included in. If None get the first index.

                **Example:**

                .. code-block:: python

                    enthalpy = engine.standardMolarEnthalpy("SiO2")
                    print(enthalpy)
                )doc")
        .def("standardMolarEnthalpy_i", &ChemicalEngine::standardMolarEnthalpy_i,
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar enthalpy of a specific species by its index (J/mol).
                Access specified element without bounds checking.

                :param int ispecies: Index of the species.

                **Example:**

                .. code-block:: python

                    enthalpy = engine.standardMolarEnthalpy_i(0)
                    print(enthalpy)
                )doc")

        .def("standardMolarVolume", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::standardMolarVolume),
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar volume of a species by its index. (m³/mol).
                Access specified element without bounds checking.

                :param int ispecies: Index of the species.

                **Example:**

                .. code-block:: python

                    volume = engine.standardMolarVolume(0)
                    print(volume)
                )doc")
        .def("standardMolarVolume", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::standardMolarVolume),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
                Returns the standard molar volume of a species by its index. (m³/mol).

                :param str species: Name of the species.
                :param str phase: Name of the phase the species was included in. If None get the first index.

                **Example:**

                .. code-block:: python

                    volume = engine.standardMolarVolume("SiO2")
                    print(volume)
                )doc")
        .def("standardMolarVolume_i", &ChemicalEngine::standardMolarVolume_i,
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar volume of a species by its index. (m³/mol).
                Access specified element without bounds checking.

                :param int ispecies: Index of the species.

                **Example:**

                .. code-block:: python

                    volume = engine.standardMolarVolume_i(0)
                    print(volume)
                )doc")

        .def("standardMolarEntropy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::standardMolarEntropy),
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar entropy of a specific species by its index (J/K/mol).
                Access specified element with bounds checking.

                :param int ispecies: Index of the species.
                
                **Example:**
                
                .. code-block:: python
                
                    entropy = engine.standardMolarEntropy(1)
                    print(entropy)
                )doc")
        .def("standardMolarEntropy", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::standardMolarEntropy),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
                Returns the standard molar entropy of a specific species by its index (J/K/mol).

                :param str species: Name of the species.
                :param str phase: Name of the phase the species was included in. If None get the first index.

                **Example:**

                .. code-block:: python

                    entropy = engine.standardMolarEntropy("SiO2")
                    print(entropy)
                )doc")
        .def("standardMolarEntropy_i", &ChemicalEngine::standardMolarEntropy_i,
             py::arg("ispecies"),
             R"doc(
                Returns the standard molar entropy of a specific species by its index (J/K/mol).
                Access specified element without bounds checking.

                :param int ispecies: Index of the species.

                **Example:**

                .. code-block:: python

                    entropy = engine.standardMolarEntropy_i(1)
                    print(entropy)
                )doc")

        //   .def("standardMolarInternalEnergy", &ChemicalEngine::standardMolarInternalEnergy)
        //   .def("standardMolarHelmholtzEnergy", &ChemicalEngine::standardMolarHelmholtzEnergy)
        .def("standardMolarHeatCapacityConstP",static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::standardMolarHeatCapacityConstP),
             py::arg("ispecies"),
             R"doc(
        Returns the standard molar heat capacity at constant pressure for species (J/K/mol).
        Access specified element with bounds checking.

        :param int ispecies: Index of the species.

        **Example:**
        
        .. code-block:: python
        
            heat_capacity = engine.standardMolarHeatCapacityConstP(0)
            print(heat_capacity)
        )doc")
        .def("standardMolarHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(std::string, std::optional<std::string>) const>(&ChemicalEngine::standardMolarHeatCapacityConstP),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
        Returns the standard molar heat capacity at constant pressure for species (J/K/mol).

       :param str species: Name of the species.
       :param str phase: Name of the phase the species was included in. If None get the first index.

        **Example:**

        .. code-block:: python

            heat_capacity = engine.standardMolarHeatCapacityConstP("SiO2")
            print(heat_capacity)
        )doc")
        .def("standardMolarHeatCapacityConstP_i", &ChemicalEngine::standardMolarHeatCapacityConstP_i,
             py::arg("ispecies"),
             R"doc(
        Returns the standard molar heat capacity at constant pressure for species (J/K/mol).
        Access specified element without bounds checking.

        :param int ispecies: Index of the species.

        **Example:**

        .. code-block:: python

            heat_capacity = engine.standardMolarHeatCapacityConstP_i(0)
            print(heat_capacity)
        )doc")
        //   .def("standardMolarHeatCapacityConstV", &ChemicalEngine::standardMolarHeatCapacityConstV)
        .def("phaseMolarGibbsEnergy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseMolarGibbsEnergy),
             py::arg("iphase"),
             R"doc(
        Returns the molar Gibbs energy of a specific phase by its index (J/mol).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_gibbs_energy = engine.phaseMolarGibbsEnergy(0)
            print(molar_gibbs_energy)
        )doc")
        .def("phaseMolarGibbsEnergy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseMolarGibbsEnergy),
             py::arg("phase"),
             R"doc(
        Returns the molar Gibbs energy of a specific phase by its index (J/mol).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            molar_gibbs_energy = engine.phaseMolarGibbsEnergy("aq_gen")
            print(molar_gibbs_energy)
        )doc")
        .def("phaseMolarGibbsEnergy_i", &ChemicalEngine::phaseMolarGibbsEnergy_i,
             py::arg("iphase"),
             R"doc(
        Returns the molar Gibbs energy of a specific phase by its index (J/mol).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_gibbs_energy = engine.phaseMolarGibbsEnergy_i(0)
            print(molar_gibbs_energy)
        )doc")

        .def("phaseMolarEnthalpy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseMolarEnthalpy),
             py::arg("iphase"),
             R"doc(
        Returns the molar enthalpy of a specific phase by its index (J/mol).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_enthalpy = engine.phaseMolarEnthalpy(0)
            print(molar_enthalpy)
        )doc")
        .def("phaseMolarEnthalpy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseMolarEnthalpy),
             py::arg("phase"),
             R"doc(
        Returns the molar enthalpy of a specific phase (J/mol).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            molar_enthalpy = engine.phaseMolarEnthalpy("aq_gen")
            print(molar_enthalpy)
        )doc")
        .def("phaseMolarEnthalpy_i", &ChemicalEngine::phaseMolarEnthalpy_i,
             py::arg("iphase"),
             R"doc(
        Returns the molar enthalpy of a specific phase by its index (J/mol).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_enthalpy = engine.phaseMolarEnthalpy_i(0)
            print(molar_enthalpy)
        )doc")

        .def("phaseMolarVolume", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseMolarVolume),
             py::arg("iphase"),
             R"doc(
        Returns the molar volume of a specific phase by its index (m³/mol).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            molar_volume = engine.phaseMolarVolume(0)
            print(molar_volume)
        )doc")
        .def("phaseMolarVolume", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseMolarVolume),
             py::arg("phase"),
             R"doc(
        Returns the molar volume of a specific phase (m³/mol).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            molar_volume = engine.phaseMolarVolume("aq_gen")
            print(molar_volume)
        )doc")
        .def("phaseMolarVolume_i", &ChemicalEngine::phaseMolarVolume_i,
             py::arg("iphase"),
             R"doc(
        Returns the molar volume of a specific phase by its index (m³/mol).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_volume = engine.phaseMolarVolume_i(0)
            print(molar_volume)
        )doc")

        .def("phaseMolarEntropy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseMolarEntropy),
             py::arg("iphase"),
             R"doc(
        Returns the molar entropy of a specific phase by its index (J/K/mol).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_entropy = engine.phaseMolarEntropy(0)
            print(molar_entropy)
        )doc")
        .def("phaseMolarEntropy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseMolarEntropy),
             py::arg("phase"),
             R"doc(
        Returns the molar entropy of a specific phase (J/K/mol).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            molar_entropy = engine.phaseMolarEntropy("aq_gen")
            print(molar_entropy)
        )doc")
        .def("phaseMolarEntropy_i", &ChemicalEngine::phaseMolarEntropy_i,
             py::arg("iphase"),
             R"doc(
        Returns the molar entropy of a specific phase by its index (J/K/mol).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            molar_entropy = engine.phaseMolarEntropy_i(0)
            print(molar_entropy)
        )doc")

        //    .def("phaseMolarInternalEnergy", &ChemicalEngine::phaseMolarInternalEnergy)
        //    .def("phaseMolarHelmholtzEnergy", &ChemicalEngine::phaseMolarHelmholtzEnergy)
        .def("phaseMolarHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseMolarHeatCapacityConstP),
             py::arg("iphase"),
             R"doc(
        Returns the molar heat capacity at constant pressure of a specific species by its index (J/K/mol).
        Access specified element with bounds checking.

        :param int iphase: Index of the species.

        **Example:**

        .. code-block:: python

            molar_heat_capacity = engine.phaseMolarHeatCapacityConstP(0)
            print(molar_heat_capacity)
        )doc")
        .def("phaseMolarHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseMolarHeatCapacityConstP),
             py::arg("phase"),
             R"doc(
        Returns the molar heat capacity at constant pressure of a specific species (J/K/mol).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            molar_heat_capacity = engine.phaseMolarHeatCapacityConstP("aq_gen")
            print(molar_heat_capacity)
        )doc")
        .def("phaseMolarHeatCapacityConstP_i", &ChemicalEngine::phaseMolarHeatCapacityConstP_i,
             py::arg("iphase"),
             R"doc(
        Returns the molar heat capacity at constant pressure of a specific species by its index (J/K/mol).
        Access specified element without bounds checking.

        :param int iphase: Index of the species.

        **Example:**

        .. code-block:: python

            molar_heat_capacity = engine.phaseMolarHeatCapacityConstP_i(0)
            print(molar_heat_capacity)
        )doc")

        //    .def("phaseMolarHeatCapacityConstV", &ChemicalEngine::phaseMolarHeatCapacitiesConstV)
        .def("phaseSpecificGibbsEnergy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseSpecificGibbsEnergy),
             py::arg("iphase"),
             R"doc(
        Returns the specific Gibbs energy of a specific phase by its index (J/kg).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_gibbs_energy = engine.phaseSpecificGibbsEnergy(0)
            print(specific_gibbs_energy)
        )doc")
        .def("phaseSpecificGibbsEnergy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseSpecificGibbsEnergy),
             py::arg("phase"),
             R"doc(
        Returns the specific Gibbs energy of a specific phase (J/kg).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            specific_gibbs_energy = engine.phaseSpecificGibbsEnergy("aq_gen")
            print(specific_gibbs_energy)
        )doc")
        .def("phaseSpecificGibbsEnergy_i", &ChemicalEngine::phaseSpecificGibbsEnergy_i,
             R"doc(
        Returns the specific Gibbs energy of a specific phase by its index (J/kg).
        Access specified element without bounds checking.

        :param int index: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_gibbs_energy = engine.phaseSpecificGibbsEnergy_i(0)
            print(specific_gibbs_energy)
        )doc")

        .def("phaseSpecificEnthalpy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseSpecificEnthalpy),
             py::arg("iphase"),
             R"doc(
        Returns the specific enthalpy of a specific phase by its index (J/kg).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            specific_enthalpy = engine.phaseSpecificEnthalpy(0)
            print(specific_enthalpy)
        )doc")
        .def("phaseSpecificEnthalpy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseSpecificEnthalpy),
             py::arg("phase"),
             R"doc(
        Returns the specific enthalpy of a specific phase (J/kg).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            specific_enthalpy = engine.phaseSpecificEnthalpy("aq_gen")
            print(specific_enthalpy)
        )doc")
        .def("phaseSpecificEnthalpy_i", &ChemicalEngine::phaseSpecificEnthalpy_i,
             py::arg("iphase"),
             R"doc(
        Returns the specific enthalpy of a specific phase by its index (J/kg).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_enthalpy = engine.phaseSpecificEnthalpy_i(0)
            print(specific_enthalpy)
        )doc")

        .def("phaseSpecificVolume", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseSpecificVolume),
             py::arg("iphase"),
             R"doc(
        Returns the specific volume of a specific phase by its index (m³/kg).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_volume = engine.phaseSpecificVolume(0)
            print(specific_volume)
        )doc")
        .def("phaseSpecificVolume", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseSpecificVolume),
             py::arg("phase"),
             R"doc(
        Returns the specific volume of a specific phase (m³/kg).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            specific_volume = engine.phaseSpecificVolume("aq_gen")
            print(specific_volume)
        )doc")
        .def("phaseSpecificVolume_i", &ChemicalEngine::phaseSpecificVolume_i,
             R"doc(
        Returns the specific volume of a specific phase by its index (m³/kg).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_volume = engine.phaseSpecificVolume_i(0)
            print(specific_volume)
        )doc")

        .def("phaseSpecificEntropy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseSpecificEntropy),
             py::arg("iphase"),
             R"doc(
        Returns the specific entropy of a specific phase by its index (J/K/kg).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            specific_entropy = engine.phaseSpecificEntropy(0)
            print(specific_entropy)
        )doc")
        .def("phaseSpecificEntropy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseSpecificEntropy),
             py::arg("phase"),
             R"doc(
        Returns the specific entropy of a specific phase (J/K/kg).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            specific_entropy = engine.phaseSpecificEntropy("aq_gen")
            print(specific_entropy)
        )doc")
        .def("phaseSpecificEntropy_i", &ChemicalEngine::phaseSpecificEntropy_i,
             py::arg("iphase"),
             R"doc(
        Returns the specific entropy of a specific phase by its index (J/K/kg).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_entropy = engine.phaseSpecificEntropy_i(0)
            print(specific_entropy)
        )doc")

        //    .def("phaseSpecificInternalEnergy", &ChemicalEngine::phaseSpecificInternalEnergy)
        //    .def("phaseSpecificHelmholtzEnergy", &ChemicalEngine::phaseSpecificHelmholtzEnergy)
        .def("phaseSpecificHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseSpecificHeatCapacityConstP),
             py::arg("iphase"),
             R"doc(
        Returns the specific heat capacity at constant pressure of a specific phase by its index (J/K/kg).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_heat_capacity = engine.phaseSpecificHeatCapacityConstP(0)
            print(specific_heat_capacity)
        )doc")
        .def("phaseSpecificHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseSpecificHeatCapacityConstP),
             py::arg("phase"),
             R"doc(
        Returns the specific heat capacity at constant pressure of a specific phase (J/K/kg).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            specific_heat_capacity = engine.phaseSpecificHeatCapacityConstP("aq_gen")
            print(specific_heat_capacity)
        )doc")
        .def("phaseSpecificHeatCapacityConstP_i", &ChemicalEngine::phaseSpecificHeatCapacityConstP_i,
             py::arg("iphase"),
             R"doc(
        Returns the specific heat capacity at constant pressure of a specific phase by its index (J/K/kg).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            specific_heat_capacity = engine.phaseSpecificHeatCapacityConstP_i(0)
            print(specific_heat_capacity)
        )doc")

        //    .def("phaseSpecificHeatCapacityConstV", &ChemicalEngine::phaseSpecificHeatCapacityConstVl)
        .def("phaseDensities", &ChemicalEngine::phaseDensities,
             R"doc(
        Returns the densities of all phases in the system (kg/m³).
        
        **Example:**
        
        .. code-block:: python
        
            densities = engine.phaseDensities()
            print(densities)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseDensity", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseDensity),
             py::arg("iphase"),
             R"doc(
        Returns the density of a specific phase by its index (kg/m³).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            density = engine.phaseDensity(0)
            print(density)
        )doc")
        .def("phaseDensity", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseDensity),
             py::arg("phase"),
             R"doc(
        Returns the density of a specific phase (kg/m³).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            density = engine.phaseDensity("aq_gen")
            print(density)
        )doc")
        .def("phaseDensity_i", &ChemicalEngine::phaseDensity_i,
             py::arg("iphase"),
             R"doc(
        Returns the density of a specific phase by its index (kg/m³).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            density = engine.phaseDensity_i(0)
            print(density)
        )doc")

        .def("phaseMasses", &ChemicalEngine::phaseMasses,
             R"doc(
        Returns the masses of all phases in the system (kg).
        
        **Example:**
        
        .. code-block:: python
        
            masses = engine.phaseMasses()
            print(masses)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseMass", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseMass),
             py::arg("iphase"),
             R"doc(
        Returns the mass of a specific phase by its index (kg).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            mass = engine.phaseMass(0)
            print(mass)
        )doc")
        .def("phaseMass", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseMass),
             py::arg("phase"),
             R"doc(
        Returns the mass of a specific phase (kg).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            mass = engine.phaseMass("aq_gen")
            print(mass)
        )doc")
        .def("phaseMass_i", &ChemicalEngine::phaseMass_i,
             py::arg("iphase"),
             R"doc(
        Returns the mass of a specific phase by its index (kg).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            mass = engine.phaseMass_i(0)
            print(mass)
        )doc")

        .def("phaseAmounts", &ChemicalEngine::phaseAmounts,
             R"doc(
        Returns the amounts of all phases in the system (mol).
        
        **Example:**
        
        .. code-block:: python
        
            amounts = engine.phaseAmounts()
            print(amounts)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseAmount", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseAmount),
             py::arg("iphase"),
             R"doc(
        Returns the amount of a specific phase by its index (mol).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            amount = engine.phaseAmount(0)
            print(amount)
        )doc")
        .def("phaseAmount", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseAmount),
             py::arg("phase"),
             R"doc(
        Returns the amount of a specific phase (mol).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            amount = engine.phaseAmount("aq_gen")
            print(amount)
        )doc")
        .def("phaseAmount_i", &ChemicalEngine::phaseAmount_i,
             py::arg("iphase"),
             R"doc(
        Returns the amount of a specific phase by its index (mol).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            amount = engine.phaseAmount_i(0)
            print(amount)
        )doc")

        .def("phaseVolumes", &ChemicalEngine::phaseVolumes,
             R"doc(
        Returns the volumes of all phases in the system (m³).
        
        **Example:**
        
        .. code-block:: python
        
            volumes = engine.phaseVolumes()
            print(volumes)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseVolume", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseVolume),
             py::arg("iphase"),
             R"doc(
        Returns the volume of a specific phase by its index (m³).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            volume = engine.phaseVolume(0)
            print(volume)
        )doc")
        .def("phaseVolume", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseVolume),
             py::arg("phase"),
             R"doc(
        Returns the volume of a specific phase (m³).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            volume = engine.phaseVolume("aq_gen")
            print(volume)
        )doc")
        .def("phaseVolume_i", &ChemicalEngine::phaseVolume_i,
             py::arg("iphase"),
             R"doc(
        Returns the iphase of a specific phase by its index (m³).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            volume = engine.phaseVolume_i(0)
            print(volume)
        )doc")

        .def("phaseEnthalpies", &ChemicalEngine::phaseEnthalpies,
             R"doc(
        Returns the enthalpies of all phases in the system (J).
        
        **Example:**
        
        .. code-block:: python
        
            enthalpies = engine.phaseEnthalpies()
            print(enthalpies)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseEnthalpy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseEnthalpy),
             py::arg("iphase"),
             R"doc(
        Returns the enthalpy of a specific phase by its index (J).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            enthalpy = engine.phaseEnthalpy(0)
            print(enthalpy)
        )doc")
        .def("phaseEnthalpy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseEnthalpy),
             py::arg("phase"),
             R"doc(
        Returns the enthalpy of a specific phase (J).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            enthalpy = engine.phaseEnthalpy("aq_gen")
            print(enthalpy)
        )doc")
        .def("phaseEnthalpy_i", &ChemicalEngine::phaseEnthalpy_i,
             py::arg("iphase"),
             R"doc(
        Returns the enthalpy of a specific phase by its index (J).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            enthalpy = engine.phaseEnthalpy_i(0)
            print(enthalpy)
        )doc")

        .def("phaseEntropies", &ChemicalEngine::phaseEntropies,
             R"doc(
        Returns the entropies of all phases in the system (J/K).
        
        **Example:**
        
        .. code-block:: python
        
            entropies = engine.phaseEntropies()
            print(entropies)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseEntropy", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseEntropy),
             py::arg("iphase"),
             R"doc(
        Returns the entropy of a specific phase by its index (J/K).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            entropy = engine.phaseEntropy(0)
            print(entropy)
        )doc")
        .def("phaseEntropy", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseEntropy),
             py::arg("phase"),
             R"doc(
        Returns the entropy of a specific phase (J/K).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            entropy = engine.phaseEntropy("aq_gen")
            print(entropy)
        )doc")
        .def("phaseEntropy_i", &ChemicalEngine::phaseEntropy_i,
             py::arg("iphase"),
             R"doc(
        Returns the entropy of a specific phase by its index (J/K).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            entropy = engine.phaseEntropy_i(0)
            print(entropy)
        )doc")

        .def("phaseHeatCapacitiesConstP", &ChemicalEngine::phaseHeatCapacitiesConstP,
             R"doc(
        Returns the isobaric heat capacities of all phases in the system (J/K).
        
        **Example:**
        
        .. code-block:: python
        
            heat_capacities = engine.phaseHeatCapacitiesConstP()
            print(heat_capacities)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseHeatCapacityConstP),
             py::arg("iphase"),
             R"doc(
        Returns the isobaric heat capacity of a specific phase by its index (J/K).
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.
        
        **Example:**
        
        .. code-block:: python
        
            heat_capacity = engine.phaseHeatCapacityConstP(0)
            print(heat_capacity)
        )doc")
        .def("phaseHeatCapacityConstP", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseHeatCapacityConstP),
             py::arg("phase"),
             R"doc(
        Returns the isobaric heat capacity of a specific phase (J/K).

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            heat_capacity = engine.phaseHeatCapacityConstP("aq_gen")
            print(heat_capacity)
        )doc")
        .def("phaseHeatCapacityConstP_i", &ChemicalEngine::phaseHeatCapacityConstP_i,
             py::arg("iphase"),
             R"doc(
        Returns the isobaric heat capacity of a specific phase by its index (J/K).
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            heat_capacity = engine.phaseHeatCapacityConstP_i(0)
            print(heat_capacity)
        )doc")

        .def("phaseSatIndices", &ChemicalEngine::phaseSatIndices,
             R"doc(
        Returns the saturation indices of all phases in the system.
        
        **Example:**
        
        .. code-block:: python
        
            sat_indices = engine.phaseSatIndices()
            print(sat_indices)
        )doc",
             py::return_value_policy::reference_internal)

        .def("phaseSatIndex", static_cast<double (ChemicalEngine::*)(Index) const>(&ChemicalEngine::phaseSatIndex),
             py::arg("iphase"),
             R"doc(
        Returns the saturation index of a specific phase by its index.
        Access specified element with bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            sat_index = engine.phaseSatIndex(0)
            print(sat_index)
        )doc")
        .def("phaseSatIndex", static_cast<double (ChemicalEngine::*)(std::string) const>(&ChemicalEngine::phaseSatIndex),
             py::arg("phase"),
             R"doc(
        Returns the saturation index of a specific phase.

        :param str phase: Name of the phase.

        **Example:**

        .. code-block:: python

            sat_index = engine.phaseSatIndex("aq_gen")
            print(sat_index)
        )doc")
        .def("phaseSatIndex_i", &ChemicalEngine::phaseSatIndex_i,
             py::arg("iphase"),
             R"doc(
        Returns the saturation index of a specific phase by its index.
        Access specified element without bounds checking.

        :param int iphase: Index of the phase.

        **Example:**

        .. code-block:: python

            sat_index = engine.phaseSatIndex_i(0)
            print(sat_index)
        )doc")

        .def("systemMass", &ChemicalEngine::systemMass,
             R"doc(
        Returns the total mass of the system (kg).
        
        **Example:**
        
        .. code-block:: python
        
            mass = engine.systemMass()
            print(mass)
        )doc")

        .def("systemVolume", &ChemicalEngine::systemVolume,
             R"doc(
        Returns the total volume of the system (m³).
        
        **Example:**
        
        .. code-block:: python
        
            volume = engine.systemVolume()
            print(volume)
        )doc")

        .def("systemGibbsEnergy", &ChemicalEngine::systemGibbsEnergy,
             R"doc(
        Returns the total Gibbs energy of the system (J).
        
        **Example:**
        
        .. code-block:: python
        
            gibbs_energy = engine.systemGibbsEnergy()
            print(gibbs_energy)
        )doc")

        .def("systemEnthalpy", &ChemicalEngine::systemEnthalpy,
             R"doc(
        Returns the total enthalpy of the system (J).
        
        **Example:**
        
        .. code-block:: python
        
            enthalpy = engine.systemEnthalpy()
            print(enthalpy)
        )doc")

        .def("systemEntropy", &ChemicalEngine::systemEntropy,
             R"doc(
        Returns the total entropy of the system (J/K).
        
        **Example:**
        
        .. code-block:: python
        
            entropy = engine.systemEntropy()
            print(entropy)
        )doc")

        .def("systemHeatCapacityConstP", &ChemicalEngine::systemHeatCapacityConstP,
             R"doc(
        Returns the total isobaric heat capacity of the system (J/K).
        
        **Example:**
        
        .. code-block:: python
        
            heat_capacity = engine.systemHeatCapacityConstP()
            print(heat_capacity)
        )doc")

        .def("ionicStrength", &ChemicalEngine::ionicStrength,
             R"doc(
        Returns the ionic strength of the system (mol/L).
        
        **Example:**
        
        .. code-block:: python
        
            ionic_strength = engine.ionicStrength()
            print(ionic_strength)
        )doc")

        .def("pH", &ChemicalEngine::pH,
             R"doc(
        Returns the pH of the system (dimensionless).
        
        **Example:**
        
        .. code-block:: python
        
            ph_value = engine.pH()
            print(ph_value)
        )doc")

        .def("pe", &ChemicalEngine::pe,
             R"doc(
        Returns the pe (electron activity) of the system (dimensionless).
        
        **Example:**
        
        .. code-block:: python
        
            pe_value = engine.pe()
            print(pe_value)
        )doc")

        .def("Eh", &ChemicalEngine::Eh,
             R"doc(
        Returns the redox potential (Eh) of the system (V).
        
        **Example:**
        
        .. code-block:: python
        
            eh_value = engine.Eh()
            print(eh_value)
        )doc")

        .def("aqueousPhaseName", &ChemicalEngine::aqueousPhaseName,
              R"doc(
        Returns the aqueous phase name. If empty, the aqueous phase is not in system.

        **Example:**

        .. code-block:: python

            name = engine.aqueousPhaseName()
            print(name)
        )doc")

        .def("gasPhaseName", &ChemicalEngine::gasPhaseName,
             R"doc(
        Returns the gaseous phase name. If empty, the gaseous phase is not in system.

        **Example:**

        .. code-block:: python

            name = engine.gasPhaseName()
            print(name)
        )doc")


        .def("__repr__", [](const ChemicalEngine &self)
             {
            std::stringstream ss;
            ss << self;
            return ss.str(); }, R"doc(
        Returns a string representation of the ChemicalEngine instance.
        
        **Example:**
        
        .. code-block:: python
        
            print(engine)
        )doc");
}
