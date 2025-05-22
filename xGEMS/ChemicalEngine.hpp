/**
 * @file ChemicalEngine.hpp
 * @brief Header file for the ChemicalEngine class.
 *
 * xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
 * This file defines the classes and functions for Gibbs-energy minimization
 * thermodynamic modeling with xGEMS. Units used throughout the API are:
 * - Temperature: Kelvin (K)
 * - Pressure: Pascals (Pa)
 * - Amounts: moles (mol)
 * - Mass: kilograms (kg)
 * - Volume: cubic meters (m³)
 * - Energies: Joules (J) or Joules per mole (J/mol) as appropriate.
 *
 * Example applications of these functions can be found in the demos at:
 * https://bitbucket.org/gems4/xgems/src/master/demos/
 *
 * @note The examples below are indicative; please ensure that your vectors
 *       and matrices (such as VectorConstRef) are properly defined (e.g., using Eigen).
 *
 * @author Allan Leal, Dmitrii Kulik, G.D. Miron
 * @date 2025
 *
 * license GNU General Public License v3 or later
 */

 #pragma once

 // C++ includes
 #include <iostream>
 #include <memory>
 #include <string>
 
 // xGEMS includes (assumed to define Index and Eigen-based vector/matrix types)
 #include <xGEMS/Index.hpp>
 #include <xGEMS/Eigen.hpp>
 
 namespace xGEMS {
 
 /**
  * @brief Updates logger settings.
  *
  * Configures logging behavior for the ChemicalEngine.
  *
  * @param use_cout If true, prints log messages to stdout.
  * @param logfile_name The filename for a rotating log file. If empty, file logging is disabled.
  * @param log_level Verbosity level for logging (trace = 0, debug = 1, info = 2, warn = 3, err = 4, critical = 5, off = 6 ).
  *
  * @code
  * // Example: Enable stdout logging and log to "chemical_log.txt" with verbosity level 2.
  * xGEMS::update_loggers(true, "chemical_log.txt", 2);
  * @endcode
  */
 void update_loggers(bool use_cout, const std::string& logfile_name, size_t log_level);
 
 /**
  * @brief Options for configuring a ChemicalEngine instance.
  *
  * This structure contains settings (e.g., initial approximation method) that affect
  * the performance and accuracy of the equilibrium computations.
  *
  * @note The member `warmstart` indicates whether a "smart" (fast but possibly less
  *       accurate) initial approximation is used.
  *
  * @code
  * // Example: Disable warm start when a more robust (cold start) solution is needed.
  * xGEMS::ChemicalEngineOptions options;
  * options.warmstart = false;
  * @endcode
  */
 struct ChemicalEngineOptions
 {
     /// Use smart start initial approximation (faster convergence, may be less accurate).
     bool warmstart = false;
 };
 
 /**
  * @class ChemicalEngine
  * @brief Class for equilibrium computations and thermodynamic analysis.
  *
  * The ChemicalEngine performs Gibbs energy minimization to compute the equilibrium
  * state of a chemical system. Its API provides methods to load system data, update
  * component amounts, set options, and query resulting thermodynamic properties.
  */
 class ChemicalEngine
 {
 public:

 // Initialization and Configuration
     /**
      * @brief Default constructor.
      *
      * Creates an empty ChemicalEngine instance. Use one of the initialization routines
      * (such as initialize() or initializeFromJsonStrings()) to configure the system.
      *
      * @code
      * // Example: Create an engine and then initialize it with a system file.
      * xGEMS::ChemicalEngine engine;
      * engine.initialize("my-system-dat.lst"); // Temperature in K, Pressure in Pa.
      * @endcode
      */
     ChemicalEngine();
 
     /**
      * @brief Constructs a ChemicalEngine instance by loading a GEM-Selektor project file.
      *
      * @param filename The file path for the chemical system definition (e.g., "my-system-dat.lst").
      *
      * @code
      * // Example: Directly initialize the engine from a project file.
      * xGEMS::ChemicalEngine engine("my-system-dat.lst");
      * @endcode
      */
     ChemicalEngine(std::string filename);
 
     /**
      * @brief Destructor.
      *
      * Cleans up resources. No explicit example is necessary.
      */
     virtual ~ChemicalEngine();
 
     /// @brief Deleted copy constructor.
     ChemicalEngine(const ChemicalEngine& other) = delete;
 
     /// @brief Deleted assignment operator.
     auto operator=(ChemicalEngine other) -> ChemicalEngine& = delete;
 
     /**
      * @brief Reallocates internal arrays.
      *
      * Should be invoked after the system composition changes to ensure proper memory
      * configuration.
      *
      * @code
      * // Example: After modifying system composition, reallocate arrays.
      * engine.reallocateEngineArrays();
      * @endcode
      */
     auto reallocateEngineArrays() -> void;
 
     /**
      * @brief Initializes the ChemicalEngine from a GEM-Selektor project file.
      *
      * Reads the file and configures the chemical system defined inside it.
      *
      * @param filename Path to the system definition file (e.g., "my-system-dat.lst").
      *
      * @code
      * // Example: Initialize using the system definition file.
      * xGEMS::ChemicalEngine engine;
      * engine.initialize("my-system-dat.lst");
      * @endcode
      */
     auto initialize(std::string filename) -> void;
     
     /**
      * @brief Initializes the ChemicalEngine using JSON strings.
      *
      * Configures the engine from three JSON strings:
      * - dch_json: Defines the chemical system.
      * - ipm_json: Contains calculation parameters and settings.
      * - dbr_json: Specifies the input node composition.
      * These strings can be obtained from the GEM-Selektor exported GEMS3K system files as JSON (*-imp.json, *-dbc.json, *-dbr.json).
      *
      * @param dch_json JSON string describing the chemical system.
      * @param ipm_json JSON string with parameter and algorithm settings.
      * @param dbr_json JSON string containing node composition details.
      *
      * @code
      * // Example: Initialize using JSON-based input.
      * std::string dch = "{}";
      * std::string ipm = "{}";
      * std::string dbr = "{}";
      * engine.initializeFromJsonStrings(dch, ipm, dbr);
      * @endcode
      */
     auto initializeFromJsonStrings(std::string dch_json, std::string ipm_json, std::string dbr_json) -> void;
 
     /**
      * @brief Reads a DBR file from disk.
      *
      * Loads system composition and parameter data from a DBR file.
      *
      * @param filename Path to the input DBR file (e.g., "*-dbr.json/.dat").
      *
      * @code
      * // Example: Read system composition from a DBR file.
      * engine.readDbrFromFile("my-system-dbr-0-0000.json");
      * @endcode
      */
     auto readDbrFromFile(std::string filename) -> void;
 
     /**
      * @brief Reads System DBR composition from a JSON string.
      *
      * Updates the system composition using JSON information.
      *
      * @param dbr_json JSON string containing DBR composition data.
      *
      * @code
      * // Example: Update composition dynamically from a JSON string.
      * engine.readDbrFromJsonString("{}");
      * @endcode
      */
     auto readDbrFromJsonString(std::string dbr_json) -> void;
 
     /**
      * @brief Writes the current DBR to a file (key-value format).
      *
      * Saves the state of the system (after equilibrium calculations) into a DBR file.
      *
      * @param filename Path to the output DBR file (e.g., "my-system-dbr-0-0001.dat").
      *
      * @code
      * // Example: Save the current state to"my-system-dbr-0-0001.dat".
      * engine.writeDbrToFile("my-system-dbr-0-0001.dat");
      * @endcode
      */
     auto writeDbrToFile(std::string filename) -> void;
 
     /**
      * @brief Returns the current DBR as a JSON string.
      *
      * @return A JSON string representing the system composition.
      *
      * @code
      * // Example: Retrieve the system state as JSON.
      * std::string json_state = engine.writeDbrToJsonString();
      * std::cout << json_state << std::endl;
      * @endcode
      */
     auto writeDbrToJsonString() -> const std::string;
 
     /**
      * @brief Returns the number of elements in the system.
      *
      * @return (Index) Total number of elements.
      *
      * @code
      * // Example: Get the number of elements.
      * Index nElements = engine.numElements();
      * std::cout << "Elements: " << nElements << std::endl;
      * @endcode
      */
     auto numElements() const -> Index;
 
     /**
      * @brief Returns the number of species in the system.
      *
      * @return (Index) Total number of chemical species.
      *
      * @code
      * // Example: Display the number of species.
      * Index speciesCount = engine.numSpecies();
      * std::cout << "Species: " << speciesCount << std::endl;
      * @endcode
      */
     auto numSpecies() const -> Index;
 
     /**
      * @brief Returns the number of phases in the system.
      *
      * @return (Index) Total number of phases.
      *
      * @code
      * // Example: Print the number of phases.
      * Index phases = engine.numPhases();
      * std::cout << "Phases: " << phases << std::endl;
      * @endcode
      */
     auto numPhases() const -> Index;
 
     /**
      * @brief Returns the number of species in a given phase.
      *
      * @param iphase (Index) Index of the phase.
      * @return (Index) Number of species in that phase.
      *
      * @code
      * // Example: For phase with index 0.
      * Index sppInPhase = engine.numSpeciesInPhase(0);
      * std::cout << "Species in phase 0: " << sppInPhase << std::endl;
      * @endcode
      */
     auto numSpeciesInPhase(Index iphase) const -> Index;
 
     /**
      * @brief Returns the name of an element.
      *
      * @param ielement (Index) Index of the element.
      * @return (std::string) Element name.
      *
      * @code
      * // Example: Get the name of the first element.
      * std::string element = engine.elementName(0);
      * std::cout << "Element 0: " << element << std::endl;
      * @endcode
      */
     auto elementName(Index ielement) const -> std::string;
 
     /**
      * @brief Returns the name of a species.
      *
      * @param ispecies (Index) Index of the species.
      * @return (std::string) Species name.
      *
      * @code
      * // Example: Retrieve the name of species 0.
      * std::string spName = engine.speciesName(0);
      * std::cout << "Species 0: " << spName << std::endl;
      * @endcode
      */
     auto speciesName(Index ispecies) const -> std::string;
 
     /**
      * @brief Returns the electrical charge of a species.
      *
      * @param ispecies (Index) Index of the species.
      * @return (double) Charge of the species.
      *
      * @code
      * // Example: Get the charge of species with index 0.
      * double charge = engine.speciesCharge(0);
      * std::cout << "Charge: " << charge << std::endl;
      * @endcode
      */
     auto speciesCharge(Index ispecies) const -> double;
 
     /**
      * @brief Returns the name of the phase.
      *
      * @param iphase (Index) Index of the phase.
      * @return (std::string) Phase name.
      *
      * @code
      * // Example: Retrieve the name of phase 0.
      * std::string phase = engine.phaseName(0);
      * std::cout << "Phase 0: " << phase << std::endl;
      * @endcode
      */
     auto phaseName(Index iphase) const -> std::string;
 
     /**
      * @brief Returns the index of an element by name.
      *
      * @param element (std::string) Element symbol or name.
      * @return (Index) Index of the element (or total elements if not found).
      *
      * @code
      * // Example: Get the index of Oxygen ("O").
      * Index idx = engine.indexElement("O");
      * std::cout << "Index of O: " << idx << std::endl;
      * @endcode
      */
     auto indexElement(std::string element) const -> Index;
 
     /**
      * @brief Returns the index of a species by name.
      *
      * @param species (std::string) Name of the species.
      * @return (Index) Index of the species (or total species if not found).
      *
      * @code
      * // Example: Get the index of water ("H2O").
      * Index idx = engine.indexSpecies("H2O");
      * std::cout << "Index of H2O: " << idx << std::endl;
      * @endcode
      */
     auto indexSpecies(std::string species) const -> Index;
 
     /**
      * @brief Returns all indices of species matching the specified name.
      *
      * @param species (std::string) Name of the species.
      * @return (VectorXi) Vector of indices.
      *
      * @code
      * // Example: Retrieve all indices where the species "H2O" appears.
      * VectorXi allIndices = engine.indexSpeciesAll("H2O");
      * @endcode
      */
     auto indexSpeciesAll(std::string species) const -> VectorXi;
 
     /**
      * @brief Returns the index of a phase by name.
      *
      * @param phase (std::string) Name of the phase.
      * @return (Index) Index of the phase.
      *
      * @code
      * // Example: Get index for the "aq_gen" phase.
      * Index aqueousIndex = engine.indexPhase("aq_gen");
      * @endcode
      */
     auto indexPhase(std::string phase) const -> Index;
 
     /**
      * @brief Returns all indices for phases that match the given name.
      *
      * @param phase (std::string) Name of the phase.
      * @return (VectorXi) Vector of indices for matching phases.
      *
      * @code
      * // Example: Retrieve all indices for phase "SiO2".
      * VectorXi phaseIndices = engine.indexPhaseAll("SiO2");
      * @endcode
      */
     auto indexPhaseAll(std::string phase) const -> VectorXi;
 
     /**
      * @brief Returns the index of the phase containing a given species.
      *
      * @param ispecies (Index) Index of the species.
      * @return (Index) Index of the phase that contains the species.
      *
      * @code
      * // Example: Find the phase of the species with index 0.
      * Index phaseIdx = engine.indexPhaseWithSpecies(0);
      * @endcode
      */
     auto indexPhaseWithSpecies(Index ispecies) const -> Index;
 
     /**
      * @brief Returns the index of the first species in a specified phase.
      *
      * @param iphase (Index) Index of the phase.
      * @return (Index) Index of the first species in that phase.
      *
      * @code
      * // Example: Get the first species index in phase 0.
      * Index firstSpec = engine.indexFirstSpeciesInPhase(0);
      * @endcode
      */
     auto indexFirstSpeciesInPhase(Index iphase) const -> Index;
 
     /**
      * @brief Returns molar masses of elements.
      *
      * @return (VectorConstRef) Vector of molar masses (kg/mol) for each element.
      *
      * @code
      * // Example: Retrieve the element molar masses vector.
      * auto molarMasses = engine.elementMolarMasses();
      * @endcode
      */
     auto elementMolarMasses() const -> VectorConstRef;
 
     /**
      * @brief Returns molar masses of species.
      *
      * @return (VectorConstRef) Vector of species molar masses (kg/mol).
      *
      * @code
      * // Example: Get the species molar masses.
      * auto spMolarMasses = engine.speciesMolarMasses();
      * @endcode
      */
     auto speciesMolarMasses() const -> VectorConstRef;
 
     /**
      * @brief Returns the formula matrix (elements x species).
      *
      * @return (MatrixConstRef) Formula matrix where rows correspond to elements and
      * columns correspond to species.
      *
      * @code
      * // Example: Retrieve the formula matrix.
      * auto fMatrix = engine.formulaMatrix();
      * @endcode
      */
     auto formulaMatrix() const -> MatrixConstRef;
 
     /**
      * @brief Sets the complete species amounts vector.
      *
      * @param n (VectorConstRef) Vector of species amounts (mol).
      *
      * @code
      * // Example: Set all species amounts to 1.0 mol.
      * Eigen::VectorXd n(engine.numSpecies());
      * n.setOnes();
      * engine.setSpeciesAmounts(n);
      * @endcode
      */
     auto setSpeciesAmounts(VectorConstRef n) -> void;
 
     /**
      * @brief Sets the amount for a species identified by name (as in the GEMS system). 
      *
      * @param name (std::string) Species name.
      * @param amount (double) New amount in mol.
      *
      * @code
      * // Example: Set the amount of "CaSO4@" to 0.01 mol.
      * engine.setSpeciesAmount("CaSO4@", 0.01);
      * @endcode
      */
     auto setSpeciesAmount(std::string name, double amount) -> void;
 
     /**
      * @brief Sets the amount for a species identified by its index.
      *
      * @param ispecies (Index) Index of the species.
      * @param amount (double) New amount in mol.
      *
      * @code
      * // Example: Set species at index 0 to 0.01 mol.
      * engine.setSpeciesAmount(0, 0.01);
      * @endcode
      */
     auto setSpeciesAmount(Index ispecies, double amount) -> void;
 
     /**
      * @brief Sets the options for the ChemicalEngine.
      *
      * @param options (const ChemicalEngineOptions&) Configuration options.
      *
      * @code
      * // Example: Disable warm start for higher accuracy.
      * xGEMS::ChemicalEngineOptions opt;
      * opt.warmstart = false;
      * engine.setOptions(opt);
      * @endcode
      */
     auto setOptions(const ChemicalEngineOptions& options) -> void;
 
     /**
      * @brief Configures the engine to use a warm (smart) start.
      *
      * Uses previous equilibrium as initial guess (faster convergence).
      *
      * @code
      * // Example: Set smart initial approximation.
      * engine.setWarmStart();
      * @endcode
      */
     auto setWarmStart() -> void;
 
     /**
      * @brief Configures the engine to use a cold start.
      *
      * Uses a simplex LP initial guess (slower but may yield more accurate results).
      *
      * @code
      * // Example: Force a cold start.
      * engine.setColdStart();
      * @endcode
      */
     auto setColdStart() -> void;
 
     /**
      * @brief Sets an upper bound for a species identified by name.
      *
      * If amount < 0, resets to default 1e6.
      *
      * @param name (std::string) Species name.
      * @param amount (double) Upper limit in mol.
      *
      * @code
      * // Example: Set the upper limit for "SiO2" to 0.1 mol.
      * engine.setSpeciesUpperLimit("SiO2", 0.1);
      * @endcode
      */
     auto setSpeciesUpperLimit(std::string name, double amount) -> void;
 
     /**
      * @brief Sets a lower bound for a species identified by name.
      *
      * If amount < 0, resets to default 0.
      *
      * @param name (std::string) Species name.
      * @param amount (double) Lower limit in mol.
      *
      * @code
      * // Example: Set the lower limit for "SiO2" to 0.05 mol.
      * engine.setSpeciesLowerLimit("SiO2", 0.05);
      * @endcode
      */
     auto setSpeciesLowerLimit(std::string name, double amount) -> void;
 
     /**
      * @brief Sets an upper bound (maximum amount allowed to form) for a species identified by its index.
      *
      * If amount < 0, resets to default 1e6.
      *
      * @param ispecies (Index) Species index.
      * @param amount (double) Upper limit in mol.
      *
      * @code
      * // Example: Set upper limit for species index 0.
      * engine.setSpeciesUpperLimit(0, 0.1);
      * @endcode
      */
     auto setSpeciesUpperLimit(Index ispecies, double amount) -> void;
 
     /**
      * @brief Sets a lower bound (minimum amount allowed to form) for a species identified by its index.
      *
      * If amount < 0, resets to default 0.
      *
      * @param ispecies (Index) Species index.
      * @param amount (double) Lower limit in mol.
      *
      * @code
      * // Example: Set lower limit for species index 0.
      * engine.setSpeciesLowerLimit(0, 0.05);
      * @endcode
      */
     auto setSpeciesLowerLimit(Index ispecies, double amount) -> void;
 
     /**
      * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
      *
      * @param name (std::string) Species name.
      * @param value (double) Standard molar Gibbs energy in J/mol.
      *
      * @code
      * // Example: Set the standard molar Gibbs energy of H2O.
      * engine.setStandardMolarGibbsEnergy("H2O", -237140); // Value in J/mol.
      * @endcode
      */
     auto setStandardMolarGibbsEnergy(std::string name, double value) -> void;
 
     /**
      * @brief Sets upper limits for all species.
      *
      * @param n (VectorConstRef) Vector of upper limits in mol.
      *
      * @code
      * // Example: Set all species' upper limits to values provided in vector "limits".
      * engine.setSpeciesUpperLimits(limits);
      * @endcode
      */
     auto setSpeciesUpperLimits(VectorConstRef n) -> void;
 
     /**
      * @brief Sets lower limits for all species.
      *
      * @param n (VectorConstRef) Vector of lower limits in mol.
      *
      * @code
      * // Example: Set all species' lower limits.
      * engine.setSpeciesLowerLimits(lowerLimits);
      * @endcode
      */
     auto setSpeciesLowerLimits(VectorConstRef n) -> void;
 
     /**
      * @brief Sets the pressure and temperature without computing equilibrium.
      *
      * @param P Pressure in Pascals (Pa).
      * @param T Temperature in Kelvin (K).
      * @return (bool) False if PT was set correctly, true if out of range.
      *
      * @code
      * // Example: Set PT to 101325 Pa and 298.15 K.
      * bool error = engine.setPT(101325, 298.15);
      * std::cout << (error ? "PT error" : "PT set correctly") << std::endl;
      * @endcode
      */
     auto setPT( double P, double T) const -> bool;  
     
     /**
      * @brief Sets the amounts of elements (vector b) without computing equilibrium.
      *
      * @param b (VectorConstRef) Vector of element amounts in mol (order based on the index of elements in the chemical system).
      *
      * @code
      * // Example: Set elements amounts of a solution with 0.001 M MgCl2 0.001 M CO2
      * // elements order: C, Ca, Cl, H, Mg, O, Zz (charge)
      * Eigen::VectorXd b(7);
      * b << 0.001, 1e-09, 0.004, 110.6837, 0.002, 55.34405, 0;  // in mol
      * engine.setB(b);
      * @endcode
      */
     auto setB( VectorConstRef b) -> void;
 
     /**
      * @brief Re-equilibrates the system with (or without) a warm start.
      *
      * @param warmstart (bool) If true, uses previous data as an initial guess.
      * @return (int) Return code of the re-equilibration.
      *
      * @code
      * // Example: Re-equilibrate with a warm start.
      * int status = engine.reequilibrate(true);
      * @endcode
      * 
      * - 0: No GEM re-calculation needed
      * - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
      * - 2: OK after GEM calculation with LPP AIA
      * - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
      * - 4: Failure (no result) in GEM calculation with LPP AIA
      * - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
      * - 6: OK after GEM calculation with SIA
      * - 7: Bad (not fully trustful) result after GEM calculation with SIA
      * - 8: Failure (no result) in GEM calculation with SIA
      * - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.
      */
     auto reequilibrate(bool warmstart) -> int;

     /**
      * @brief Re-equilibrates the system with (or without) a warm start.
      *
      * @return (int) Return code of the re-equilibration.
      *
      * @code
      * // Example: Re-equilibrate.
      * int status = engine.reequilibrate();
      * @endcode
      * 
      * - 0: No GEM re-calculation needed
      * - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
      * - 2: OK after GEM calculation with LPP AIA
      * - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
      * - 4: Failure (no result) in GEM calculation with LPP AIA
      * - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
      * - 6: OK after GEM calculation with SIA
      * - 7: Bad (not fully trustful) result after GEM calculation with SIA
      * - 8: Failure (no result) in GEM calculation with SIA
      * - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.
      */
     auto reequilibrate() -> int;
 
     /**
      * @brief Computes the equilibrium state.
      *
      * Uses temperature (K), pressure (Pa), and element amounts (in mol) to compute equilibrium.
      *
      * @param T Temperature in Kelvin (K).
      * @param P Pressure in Pascals (Pa).
      * @param b Vector of element amounts (mol).
      * @return (int) Return code of the equilibrium solver.
      *
      * @code
      * // Example: This vector need to be defined and assembled by the user with the number of elements in the system.
      * // The order of the elements in the vector is based on the index of elements in the chemical system.
      * // Example: Set elements amounts of a solution with 0.001 M MgCl2 0.001 M CO2
      * // elements order: C, Ca, Cl, H, Mg, O, Zz (charge)
      * Eigen::VectorXd b(7);
      * b << 0.001, 1e-09, 0.004, 110.6837, 0.002, 55.34405, 0;  // in mol
      * 
      * int retcode = engine.equilibrate(298.15, 101325, b);
      * @endcode
      * - 0: No GEM re-calculation needed
      * - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
      * - 2: OK after GEM calculation with LPP AIA
      * - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
      * - 4: Failure (no result) in GEM calculation with LPP AIA
      * - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
      * - 6: OK after GEM calculation with SIA
      * - 7: Bad (not fully trustful) result after GEM calculation with SIA
      * - 8: Failure (no result) in GEM calculation with SIA
      * - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.
      */
     auto equilibrate(double T, double P, VectorConstRef b) -> int;
 
     /**
      * @brief Checks if the equilibrium calculation converged.
      *
      * @return (bool) True if converged, false otherwise.
      *
      * @code
      * // Example: Check convergence status.
      * bool isConverged = engine.converged();
      * @endcode
      */
     auto converged() const -> bool;
 
     /**
      * @brief Returns the number of iterations performed during equilibrium.
      *
      * @return (Index) The number of iterations.
      *
      * @code
      * // Example: Display the number of iterations.
      * Index iters = engine.numIterations();
      * std::cout << "Iterations: " << iters << std::endl;
      * @endcode
      */
     auto numIterations() const -> Index;
 
     /**
      * @brief Returns the elapsed time of the equilibrium computation.
      *
      * @return (double) Elapsed time in seconds.
      *
      * @code
      * // Example: Get the computation time.
      * double timeSec = engine.elapsedTime();
      * std::cout << "Time: " << timeSec << " s" << std::endl;
      * @endcode
      */
     auto elapsedTime() const -> double;
 
     /**
      * @brief Returns the current temperature of the system.
      *
      * @return (double) Temperature in Kelvin (K).
      *
      * @code
      * // Example: Get the system temperature.
      * double temp = engine.temperature();
      * std::cout << "Temperature: " << temp << " K" << std::endl;
      * @endcode
      */
     auto temperature() const -> double;
 
     /**
      * @brief Returns the current pressure of the system.
      *
      * @return (double) Pressure in Pascals (Pa).
      *
      * @code
      * // Example: Retrieve the system pressure.
      * double press = engine.pressure();
      * std::cout << "Pressure: " << press << " Pa" << std::endl;
      * @endcode
      */
     auto pressure() const -> double;
 
     /**
      * @brief Returns the amounts of the elements.
      * 
      * The order of the elements in the vector is based on the index of elements in the chemical system.
      *
      * @return (VectorConstRef) Vector of element amounts in mol.
      *
      * @code
      * // Example: Get element amounts.
      * auto eAmounts = engine.elementAmounts();
      * @endcode
      */
     auto elementAmounts() const -> VectorConstRef;

    /**
      * @brief Returns the amounts of elements in a specified phase.
      *
      * The order of the elements in the vector is based on the index of elements in the chemical system.
      *
      * @param iphase Index of the phase.
      * @return A vector containing the amounts of each element in the phase (in mol).
      *
      * **Example usage:**
      * @code
      * // Retrieve the element amounts in phase index 0.
      * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase(0);
      * std::cout << "Element amounts in phase 0: " << elementAmounts.transpose() << std::endl;
      * @endcode
      */
    auto elementAmountsInPhase(Index iphase) const -> Vector;

    /**
      * @brief Returns the amounts of elements in a specified group of species.
      *
      * The order of the elements in the vector is based on the index of elements in the chemical system.
      *
      * @param ispecies A vector containing the indices of the species whose elemental amounts should be retrieved.
      * @return A vector containing the amounts of each element corresponding to the selected species (in mol).
      *
      * **Example usage:**
      * @code
      * // Retrieve element amounts contributed by species indices {0, 2, 5}.
      * Eigen::VectorXi speciesIndices(3);
      * speciesIndices << 0, 2, 5;
      * Eigen::VectorXd elementAmounts = engine.elementAmountsInSpecies(speciesIndices);
      * std::cout << "Element amounts from selected species: " << elementAmounts.transpose() << std::endl;
      * @endcode
      */
    auto elementAmountsInSpecies(VectorXiConstRef ispecies) const -> Vector;
 
     /**
      * @brief Returns the amounts of the species.
      *
      * @return (VectorConstRef) Vector of species amounts in mol.
      *
      * @code
      * // Example: Retrieve the species amounts.
      * auto spAmounts = engine.speciesAmounts();
      * @endcode
      */
     auto speciesAmounts() const -> VectorConstRef;
 
     /**
      * @brief Returns the amount of a species by index.
      *
      * @param ispecies (Index) Index of the species.
      * @return (double) Amount in mol.
      *
      * @code
      * // Example: Get the species amount at index 0.
      * double amount = engine.speciesAmount(0);
      * @endcode
      */
     auto speciesAmount(Index ispecies) const -> double;
 
     /**
      * @brief Returns the amount of a species identified by name.
      *
      * @param name (std::string) Species name.
      * @return (double) Amount in mol.
      *
      * @code
      * // Example: Get the amount of "OH-".
      * double waterAmt = engine.speciesAmount("OH-");
      * @endcode
      */
     auto speciesAmount(std::string name) const -> double;
 
     /**
      * @brief Returns the upper limits for all species.
      *
      * @return (VectorConstRef) Vector of upper limits in mol.
      *
      * @code
      * // Example: Fetch all species upper limits.
      * auto upperLimits = engine.speciesUpperLimits();
      * @endcode
      */
     auto speciesUpperLimits() const -> VectorConstRef;
 
     /**
      * @brief Returns the lower limits for all species.
      *
      * @return (VectorConstRef) Vector of lower limits in mol.
      *
      * @code
      * // Example: Get the species lower limits.
      * auto lowerLimits = engine.speciesLowerLimits();
      * @endcode
      */
     auto speciesLowerLimits() const -> VectorConstRef;
 
     /**
      * @brief Returns the molalities of the species.
      *
      * Assumes the aqueous phase is the first phase.
      *
      * @return (VectorConstRef) Species molalities (mol/kg H2O@).
      *
      * @code
      * // Example: Retrieve species molalities.
      * auto molalities = engine.speciesMolalities();
      * @endcode
      */
     auto speciesMolalities() const -> VectorConstRef;
 
     /**
      * @brief Returns the mole fractions of the species.
      *
      * @return (VectorConstRef) Mole fractions.
      *
      * @code
      * // Example: Get species mole fractions.
      * auto moleFrac = engine.moleFractions();
      * @endcode
      */
     auto moleFractions() const -> VectorConstRef;
 
     /**
      * @brief Returns the ln activity coefficients of the species (mole fraction scale).
      *
      * @return (VectorConstRef) ln Activity coefficients.
      *
      * @code
      * // Example: Retrieve ln activity coefficients.
      * auto lnActCoeff = engine.lnActivityCoefficients();
      * @endcode
      */
     auto lnActivityCoefficients() const -> VectorConstRef;
 
     /**
      * @brief Returns the ln activities of the species.
      *
      * @return (VectorConstRef) ln Activities.
      *
      * @code
      * // Example: Fetch ln activities.
      * auto lnActivities = engine.lnActivities();
      * @endcode
      */
     auto lnActivities() const -> VectorConstRef;
 
     /**
      * @brief Returns the ln of species concentrations.
      *
      * @return (VectorConstRef) ln concentrations.
      *
      * @code
      * // Example: Get ln concentrations.
      * auto lnConcs = engine.lnConcentrations();
      * @endcode
      */
     auto lnConcentrations() const -> VectorConstRef;
 
     /**
      * @brief Returns the chemical potentials of the species.
      *
      * @return (VectorConstRef) Chemical potentials in J/mol.
      *
      * @code
      * // Example: Retrieve species chemical potentials.
      * auto chemPot = engine.chemicalPotentials();
      * @endcode
      */
     auto chemicalPotentials() const -> VectorConstRef;
 
     /**
      * @brief Returns the standard molar Gibbs energy of a species.
      *
      * @param ispecies (Index) Species index.
      * @return (double) Standard molar Gibbs energy in J/mol.
      *
      * @code
      * // Example: For species 0.
      * double stdGibbs = engine.standardMolarGibbsEnergy(0);
      * @endcode
      */
     auto standardMolarGibbsEnergy(Index ispecies) const -> double;
 
     /**
      * @brief Returns the standard molar enthalpy of a species.
      *
      * @param ispecies (Index) Species index.
      * @return (double) Standard molar enthalpy in J/mol.
      *
      * @code
      * double stdEnthalpy = engine.standardMolarEnthalpy(0);
      * @endcode
      */
     auto standardMolarEnthalpy(Index ispecies) const -> double;
 
     /**
      * @brief Returns the standard molar volume of a species.
      *
      * @param ispecies (Index) Species index.
      * @return (double) Standard molar volume in m³/mol.
      *
      * @code
      * double stdVolume = engine.standardMolarVolume(0);
      * @endcode
      */
     auto standardMolarVolume(Index ispecies) const -> double;
     
     /**
      * @brief Returns the standard molar entropy of a species.
      *
      * @param ispecies (Index) Species index.
      * @return (double) Standard molar entropy in J/(mol·K).
      *
      * @code
      * double stdEntropy = engine.standardMolarEntropy(0);
      * @endcode
      */
     auto standardMolarEntropy(Index ispecies) const -> double;
 
     /**
      * @brief Returns the standard molar isobaric heat capacity of a species.
      *
      * @param ispecies (Index) Species index.
      * @return (double) Heat capacity in J/(mol·K).
      *
      * @code
      * double heatCap = engine.standardMolarHeatCapacityConstP(0);
      * @endcode
      */
     auto standardMolarHeatCapacityConstP(Index ispecies) const -> double;
 
     /**
      * @brief Returns the molar Gibbs energy of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Molar Gibbs energy in J/mol.
      *
      * @code
      * double phaseGibbs = engine.phaseMolarGibbsEnergy(0);
      * @endcode
      */
     auto phaseMolarGibbsEnergy(Index iphase) const -> double;
 
     /**
      * @brief Returns the molar enthalpy of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Molar enthalpy in J/mol.
      *
      * @code
      * double phaseEnthalpy = engine.phaseMolarEnthalpy(0);
      * @endcode
      */
     auto phaseMolarEnthalpy(Index iphase) const -> double;
 
     /**
      * @brief Returns the molar volume of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Molar volume in m³/mol.
      *
      * @code
      * double phaseVol = engine.phaseMolarVolume(0);
      * @endcode
      */
     auto phaseMolarVolume(Index iphase) const -> double;
 
     /**
      * @brief Returns the molar entropy of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Molar entropy in J/(mol·K).
      *
      * @code
      * double phaseEnt = engine.phaseMolarEntropy(0);
      * @endcode
      */
     auto phaseMolarEntropy(Index iphase) const -> double;
 
     /**
      * @brief Returns the molar isobaric heat capacity of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Heat capacity in J/(mol·K).
      *
      * @code
      * double phaseCp = engine.phaseMolarHeatCapacityConstP(0);
      * @endcode
      */
     auto phaseMolarHeatCapacityConstP(Index iphase) const -> double;
 
     /**
      * @brief Returns the specific Gibbs energy of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Specific Gibbs energy in J/kg.
      *
      * @code
      * double specGibbs = engine.phaseSpecificGibbsEnergy(0);
      * @endcode
      */
     auto phaseSpecificGibbsEnergy(Index iphase) const -> double;
 
     /**
      * @brief Returns the specific enthalpy of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Specific enthalpy in J/kg.
      *
      * @code
      * double specEnthalpy = engine.phaseSpecificEnthalpy(0);
      * @endcode
      */
     auto phaseSpecificEnthalpy(Index iphase) const -> double;
 
     /**
      * @brief Returns the specific volume of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Specific volume in m³/kg.
      *
      * @code
      * double specVolume = engine.phaseSpecificVolume(0);
      * @endcode
      */
     auto phaseSpecificVolume(Index iphase) const -> double;
     
     /**
      * @brief Returns the specific entropy of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Specific entropy in J/(kg·K).
      *
      * @code
      * double specEntropy = engine.phaseSpecificEntropy(0);
      * @endcode
      */
     auto phaseSpecificEntropy(Index iphase) const -> double;
 
     /**
      * @brief Returns the specific isobaric heat capacity of a phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Heat capacity in J/(kg·K).
      *
      * @code
      * double specCp = engine.phaseSpecificHeatCapacityConstP(0);
      * @endcode
      */
     auto phaseSpecificHeatCapacityConstP(Index iphase) const -> double;
 
     /**
      * @brief Returns the densities of all phases.
      *
      * @return (VectorConstRef) Vector of densities in kg/m³.
      *
      * @code
      * auto densities = engine.phaseDensities();
      * @endcode
      */
     auto phaseDensities() const -> VectorConstRef;
 
     /**
      * @brief Returns the density of a specified phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Density in kg/m³.
      *
      * @code
      * double density = engine.phaseDensity(0);
      * @endcode
      */
     auto phaseDensity(Index iphase) const -> double;
 
     /**
      * @brief Returns the masses of all phases.
      *
      * @return (VectorConstRef) Vector of masses in kg.
      *
      * @code
      * auto masses = engine.phaseMasses();
      * @endcode
      */
     auto phaseMasses() const -> VectorConstRef;
  
     /**
      * @brief Returns the mass of a specified phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Mass in kg.
      *
      * @code
      * double mass = engine.phaseMass(0);
      * @endcode
      */
     auto phaseMass(Index iphase) const -> double;
     
     /**
      * @brief Returns the molar amounts of all phases.
      *
      * @return (VectorConstRef) Vector of phase amounts in mol.
      *
      * @code
      * auto phaseMol = engine.phaseAmounts();
      * @endcode
      */
     auto phaseAmounts() const -> VectorConstRef;
 
     /**
      * @brief Returns the molar amount of a specific phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Amount in mol.
      *
      * @code
      * double amt = engine.phaseAmount(0);
      * @endcode
      */
     auto phaseAmount(Index iphase) const -> double;
 
     /**
      * @brief Returns the volumes of all phases.
      *
      * @return (VectorConstRef) Vector of volumes in m³.
      *
      * @code
      * auto volumes = engine.phaseVolumes();
      * @endcode
      */
     auto phaseVolumes() const -> VectorConstRef;
  
     /**
      * @brief Returns the volume of a specific phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Volume in m³.
      *
      * @code
      * double vol = engine.phaseVolume(0);
      * @endcode
      */
     auto phaseVolume(Index iphase) const -> double;  
 
     /**
      * @brief Returns the enthalpies of all phases.
      *
      * @return (VectorConstRef) Vector of enthalpies in J.
      *
      * @code
      * auto enth = engine.phaseEnthalpies();
      * @endcode
      */
     auto phaseEnthalpies() const -> VectorConstRef;
     
     /**
      * @brief Returns the enthalpy of a specific phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Enthalpy in J.
      *
      * @code
      * double enth0 = engine.phaseEnthalpy(0);
      * @endcode
      */
     auto phaseEnthalpy(Index iphase) const -> double;
 
     /**
      * @brief Returns the entropies of all phases.
      *
      * @return (VectorConstRef) Vector of entropies in J/K.
      *
      * @code
      * auto ents = engine.phaseEntropies();
      * @endcode
      */
     auto phaseEntropies() const -> VectorConstRef;
     
     /**
      * @brief Returns the entropy of a specified phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Entropy in J/(K).
      *
      * @code
      * double ent0 = engine.phaseEntropy(0);
      * @endcode
      */
     auto phaseEntropy(Index iphase) const -> double;
 
     /**
      * @brief Returns the heat capacities (Cp) of all phases.
      *
      * @return (VectorConstRef) Vector of heat capacities in J/K.
      *
      * @code
      * auto cps = engine.phaseHeatCapacitiesConstP();
      * @endcode
      */
     auto phaseHeatCapacitiesConstP() const -> VectorConstRef;
     
     /**
      * @brief Returns the heat capacity of a specified phase.
      *
      * @param iphase (Index) Phase index.
      * @return (double) Heat capacity in J/K.
      *
      * @code
      * double cp0 = engine.phaseHeatCapacityConstP(0);
      * @endcode
      */
     auto phaseHeatCapacityConstP(Index iphase) const -> double;
 
     /**
      * @brief Returns the saturation indices of all phases (log₁₀ units).
      *
      * @return (VectorConstRef) Vector of saturation indices.
      *
      * @code
      * auto satIdx = engine.phaseSatIndices();
      * @endcode
      */
     auto phaseSatIndices() const -> VectorConstRef;
        
     /**
      * @brief Returns the saturation index of a specific phase (log₁₀ units).
      *
      * @param iphase (Index) Phase index.
      * @return (double) Saturation index.
      *
      * @code
      * double sat0 = engine.phaseSatIndex(0);
      * @endcode
      */
     auto phaseSatIndex(Index iphase) const -> double;
 
     /**
      * @brief Returns the total mass of the system.
      *
      * @return (double) System mass in kg.
      *
      * @code
      * double sysMass = engine.systemMass();
      * std::cout << "System mass: " << sysMass << " kg" << std::endl;
      * @endcode
      */
     auto systemMass() const -> double;
 
     /**
      * @brief Returns the total volume of the system.
      *
      * @return (double) Volume in m³.
      *
      * @code
      * double sysVol = engine.systemVolume();
      * std::cout << "System volume: " << sysVol << " m³" << std::endl;
      * @endcode
      */
     auto systemVolume() const -> double;
 
     /**
      * @brief Returns the ionic strength of the aqueous phase.
      *
      * @return (double) Ionic strength in molal.
      *
      * @code
      * double ionicStr = engine.ionicStrength();
      * @endcode
      */
     auto ionicStrength() const -> double;
 
     /**
      * @brief Returns the pH of the aqueous phase.
      *
      * @return (double) pH (in the activity scale (-log10 molal)).
      *
      * @code
      * double ph = engine.pH();
      * @endcode
      */
     auto pH() const -> double;
 
     /**
      * @brief Returns the pe of the aqueous phase.
      *
      * @return (double) pe (in the activity scale (-log10 molal)).
      *
      * @code
      * double pe = engine.pe();
      * @endcode
      */
     auto pe() const -> double;
 
     /**
      * @brief Returns the Eh of the aqueous phase.
      *
      * @return (double) Eh (V).
      *
      * @code
      * double eh = engine.Eh();
      * @endcode
      */
     auto Eh() const -> double;
 
     /**
      * @brief Returns the total Gibbs energy of the system.
      *
      * @return (double) Gibbs energy in J/mol.
      *
      * @code
      * double sysGibbs = engine.systemGibbsEnergy();
      * @endcode
      */
     auto systemGibbsEnergy() const -> double;
 
     /**
      * @brief Returns the total enthalpy of the system.
      *
      * @return (double) Enthalpy in J.
      *
      * @code
      * double sysEnthalpy = engine.systemEnthalpy();
      * @endcode
      */
     auto systemEnthalpy() const -> double;
 
     /**
      * @brief Returns the total entropy of the system.
      *
      * @return (double) Entropy in J/K.
      *
      * @code
      * double sysEntropy = engine.systemEntropy();
      * @endcode
      */
     auto systemEntropy() const -> double;
 
     /**
      * @brief Returns the total isobaric heat capacity of the system.
      *
      * @return (double) Heat capacity in J/K.
      *
      * @code
      * double sysCp = engine.systemHeatCapacityConstP();
      * @endcode
      */
     auto systemHeatCapacityConstP() const -> double;
 
 private:
     struct Impl;
     std::unique_ptr<Impl> pimpl;
 };
 
 /**
  * @brief Outputs the state of the ChemicalEngine.
  *
  * Overloaded stream operator for convenient printing of engine details.
  *
  * @param out Output stream.
  * @param engine ChemicalEngine instance.
  * @return (std::ostream&) Updated output stream.
  *
  * @code
  * // Example: Print the engine state.
  * std::cout << engine << std::endl;
  * @endcode
  */
 auto operator<<(std::ostream& out, const ChemicalEngine& engine) -> std::ostream&;
 
 } // namespace xGEMS