/**
 * @file ChemicalEngineMaps.hpp
 * @brief Header file for the ChemicalEngineMaps class.
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
 * /demos/
 *
 *
 * @author R.A.Patel, Dmtrii Kulik, G.D. Miron, S.Dmytriieva
 * @date 2025
 *
 * license GNU General Public License v3 or later
 */

#pragma once

#include <vector>
#include <map>
#include "ChemicalEngine.hpp"

namespace xGEMS {

/**
 * ValuesMap is a sorted associative dictionary that contains name-value pairs.
 */
using ValuesMap = std::map<std::string, double>;
/**
 * PhaseValuesMap is a sorted associative container that contains dictionaries for all phase species.
 */
using PhaseValuesMap = std::map<std::string, ValuesMap>;

/**
   * @class ChemicalEngineMaps
   * @brief Class for equilibrium computations and thermodynamic analysis using dictionaries.
   *
   * The ChemicalEngineMaps is a more convenient wrapper for Gibbs energy minimization to compute
   * the equilibrium state of a chemical system with a Pythonic naming convention and dictionaries.
   * Its API provides methods to load system data, update component amounts, and query resulting
   * thermodynamic properties as a dictionary map.
   * Gems interface in calculator format for easy  using dictionaries.
*/
class ChemicalEngineMaps
{
public:

    /**
     * @brief Constructs a ChemicalEngineMaps instance by loading a GEM-Selektor project file.
     *
     * @param input_file The file path for the chemical system definition (e.g., "my-system-dat.lst").
     * @param reset_calc (bool) If true, clear the amounts of all elements, default false.
     * @param cold_start (bool) If true, configures the engine to use a cold start, default true.
     *
     * @code
     * xGEMS::ChemicalEngineMaps engine("my-system-dat.lst");
     * @endcode
     */
    ChemicalEngineMaps(const std::string& input_file, bool reset_calc=false, bool cold_start=true);

    /**
     * @brief Constructs a ChemicalEngineMaps instance from three JSON strings.
     *
     * Accepts the same JSON content as ChemicalEngine::initializeFromJsonStrings:
     * - dch_json: chemical system definition (*-dch.json)
     * - ipm_json: IPM parameters and settings (*-ipm.json)
     * - dbr_json: node bulk composition (*-dbr.json)
     *
     * @param dch_json JSON string for the chemical system definition.
     * @param ipm_json JSON string for IPM parameters.
     * @param dbr_json JSON string for node composition.
     * @param reset_calc (bool) If true, clear all element amounts, default false.
     * @param cold_start (bool) If true, use cold start, default true.
     *
     * @code
     * xGEMS::ChemicalEngineMaps engine(dch_json, ipm_json, dbr_json);
     * @endcode
     */
    ChemicalEngineMaps(const std::string& dch_json, const std::string& ipm_json, const std::string& dbr_json,
                       bool reset_calc=false, bool cold_start=true);

    /**
     * @brief Computes the equilibrium stateof the current system.
     *
     * Uses current temperature (K), pressure (Pa), and element amounts (in mol) to compute equilibrium.
     *
     * @return (std::string) Return result string of the equilibrium solver.
     *
     * @code
     * engine.T = 298.15;
     * engine.P = 100000.0;
     * xGEMS::ValuesMap  bulk_composition = { {"C", 1e-08}, {"Ca", 1e-08}, {"Cl", 0.002},
     *                                    {"H", 111.016746657646}, {"Mg", 0.001}, {"O", 55.5083933588231},
     *                                    {"Sn", 130.841288437146}, {"Zz", 0.0} };
     * engine.set_bulk_composition(bulk_composition);
     * std::string retcode = engine.equilibrate();
     * @endcode
     *
     * - No GEM re-calculation needed
     * - Need GEM calculation with LPP (automatic) initial approximation (AIA)
     * - OK after GEM calculation with LPP AIA
     * - Bad (not fully trustful) result after GEM calculation with LPP AIA
     * - Failure (no result) in GEM calculation with LPP AIA
     * - Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
     * - OK after GEM calculation with SIA
     * - Bad (not fully trustful) result after GEM calculation with SIA
     * - Failure (no result) in GEM calculation with SIA
     * - Terminal error in GEMS3K (e.g., memory corruption). Restart required.
     */
    auto equilibrate() -> std::string;

    /**
     * @brief Computes the equilibrium state with explicit temperature, pressure and bulk composition.
     *
     * Mirrors the ChemicalEngine::equilibrate(T, P, b) call but accepts a dictionary for the
     * bulk composition instead of a vector. Sets internal T, P, and bulk composition before computing.
     *
     * @param T_new Temperature in Kelvin (K).
     * @param P_new Pressure in Pascals (Pa).
     * @param b_dict (ValuesMap) Dictionary of element amounts in mol.
     * @param min_amount (double) Minimum amount for elements not in b_dict, default 1e-15.
     * @return (std::string) Status string of the equilibrium solver.
     *
     * @code
     * xGEMS::ValuesMap b = { {"C", 0.001}, {"H", 110.0}, {"O", 55.0}, {"Zz", 0.0} };
     * std::string retcode = engine.equilibrate(298.15, 101325.0, b);
     * @endcode
     */
    auto equilibrate(double T_new, double P_new, ValuesMap b_dict, double min_amount=1e-15) -> std::string;

    /**
     * @brief Re-equilibrates the system using the current internal state (T, P, b).
     *
     * Calls the GEM solver without resetting T, P or bulk composition. Uses the warm/cold start
     * setting currently configured on the engine.
     *
     * @return (std::string) Status string of the equilibrium solver.
     *
     * @code
     * std::string retcode = engine.reequilibrate();
     * @endcode
     */
    auto reequilibrate() -> std::string;

    /**
     * @brief Re-equilibrates the system with explicit warm/cold start control.
     *
     * @param warmstart (bool) If true, uses the previous speciation as initial guess (SIA).
     * @return (std::string) Status string of the equilibrium solver.
     *
     * @code
     * std::string retcode = engine.reequilibrate(true);
     * @endcode
     */
    auto reequilibrate(bool warmstart) -> std::string;

    /**
     * @brief Sets the pressure and temperature without computing equilibrium.
     *
     * Mirrors ChemicalEngine::setPT. Also updates the T and P member variables.
     *
     * @param P Pressure in Pascals (Pa).
     * @param T Temperature in Kelvin (K).
     * @return (bool) true if PT was set correctly, false if out of range.
     *
     * @code
     * bool ok = engine.setPT(101325.0, 298.15);
     * @endcode
     */
    auto setPT(double P_new, double T_new) -> bool;

    /**
     * @brief Sets the amounts of elements (vector b) from a dictionary.
     *
     * Mirrors ChemicalEngine::setB but accepts a dictionary instead of a vector.
     * Elements absent from b_dict are left unchanged (or set to min_amount if below it).
     *
     * @param b_input (ValuesMap) Dictionary of element amounts in mol.
     * @param min_amount (double) Minimum amount for elements not provided, default 1e-15.
     *
     * @code
     * engine.setB({ {"Ca", 0.001}, {"Cl", 0.002}, {"H", 110.0}, {"O", 55.0}, {"Zz", 0.0} });
     * @endcode
     */
    auto setB(ValuesMap b_input, double min_amount=1e-15) -> void;

    /**
     * @brief Returns the current temperature of the system after the last equilibration.
     *
     * @return (double) Temperature in Kelvin (K).
     *
     * @code
     * double T = engine.temperature();
     * @endcode
     */
    auto temperature() -> double;

    /**
     * @brief Returns the current pressure of the system after the last equilibration.
     *
     * @return (double) Pressure in Pascals (Pa).
     *
     * @code
     * double P = engine.pressure();
     * @endcode
     */
    auto pressure() -> double;

    /**
     * @brief Reads a DBR file from disk, updating the system composition.
     *
     * @param filename Path to the DBR file (e.g., "*-dbr.json/.dat").
     *
     * @code
     * engine.readDbrFromFile("system-dbr-0-0000.json");
     * @endcode
     */
    auto readDbrFromFile(std::string filename) -> void;

    /**
     * @brief Reads the system DBR composition from a JSON string.
     *
     * @param dbr_json JSON string containing DBR composition data.
     *
     * @code
     * engine.readDbrFromJsonString(dbr_json_str);
     * @endcode
     */
    auto readDbrFromJsonString(std::string dbr_json) -> void;

    /**
     * @brief Writes the current DBR to a file (key-value format).
     *
     * @param filename Path to the output DBR file.
     *
     * @code
     * engine.writeDbrToFile("system-dbr-0-0001.dat");
     * @endcode
     */
    auto writeDbrToFile(std::string filename) -> void;

    /**
     * @brief Returns the current DBR as a JSON string.
     *
     * @return (std::string) JSON representation of the current system state.
     *
     * @code
     * std::string json = engine.writeDbrToJsonString();
     * @endcode
     */
    auto writeDbrToJsonString() -> const std::string;

    /**
     * @brief Configures the engine to use a cold start.
     *
     * Uses a simplex LP initial guess (slower but may yield more accurate results).
     *
     * @code
     * engine.cold_start();
     * @endcode
     */
    auto cold_start() -> void
    {
        gem.setColdStart();
    }

    /**
     * @brief Configures the engine to use a warm (smart) start.
     *
     * Uses previous equilibrium as initial guess (faster convergence).
     *
     * @code
     * engine.warm_start();
     * @endcode
     */
    void warm_start()
    {
        gem.setWarmStart();
    }

    /**
     * @brief Clear the amounts of elements (set the default amount for all components).
     *
     * @param min_amount (double) The minimum amount of element in mole, default 1e-15.
     *
     * @code
     * engine.clear(1e-15);
     * @endcode
     */
    auto clear(double min_amount=1e-15) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     *
     * @param symbol (std::string) Species symbol.
     * @param value (double) Standard molar Gibbs energy in J/mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * engine.set_species_G0("H2O", -237140); // Value in J/mol.
     * @endcode
     */
    auto set_species_G0(std::string symbol, double value, std::optional<std::string> phase= std::nullopt) -> void;

    /**
     * @brief Sets the amounts of elements (vector b).
     *
     * @param b_input (ValuesMap) Dictionary of elements amounts in mol.
     * @param min_amount (double) min amount of element in mol, default 1e-15.
     *
     * @code
     * xGEMS::ValuesMap  bulk_composition = { {"C", 1e-08}, {"Ca", 1e-08}, {"Cl", 0.002},
     *                                    {"H", 111.016746657646}, {"Mg", 0.001}, {"O", 55.5083933588231},
     *                                    {"Sn", 130.841288437146}, {"Zz", 0.0} };
     * engine.set_bulk_composition(bulk_composition);
     * @endcode
     */
    auto set_bulk_composition(ValuesMap b_input, double min_amount=1e-15) -> void;

    /**
     * @brief Removes bulk elemental aqueous solution composition from vector b.
     * Be careful as this will also remove water i.e H+ and OH-.
     *
     * @param min_amount (double) min amount of element in mol, default 1e-12.
     *
     * @code
     * engine.reset_aq_composition();
     * @endcode
     */
    auto reset_aq_composition(double min_amount=1e-12) -> void;

    /**
     * @brief Sets and gets the temperature without computing equilibrium.
     * T Temperature in Kelvin (K).
     *
     * @code
     * engine.T = 298.15;
     * @endcode
     */
    double T; // K

    /**
     * @brief Sets and gets the pressure without computing equilibrium.
     *
     * P Pressure in Pascals (Pa).
     *
     * @code
     * engine.P = 100000.0;
     * @endcode
     */
    double P; // Pa

    /**
     * @brief Returns the number of elements in the system.
     *
     * @return (Index) Total number of elements.
     *
     * @code
     * Index count = engine.nelements();
     * std::cout << "Elements: " << count << std::endl;
     * @endcode
     */
    auto nelements() -> int
    {
        return gem.numElements();
    }

    /**
     * @brief Returns the number of phases in the system.
     *
     * @return (Index) Total number of phases.
     *
     * @code
     * Index count = engine.nphases();
     * std::cout << "Phases: " << count << std::endl;
     * @endcode
     */
    auto nphases() -> int
    {
        return gem.numPhases();
    }

    /**
     * @brief Returns the number of species in the system.
     *
     * @return (Index) Total number of chemical species.
     *
     * @code
     * Index count = engine.nspecies();
     * std::cout << "Species: " << count << std::endl;
     * @endcode
     */
    auto nspecies() -> int
    {
        return gem.numSpecies();
    }

    /**
     * @brief Returns the names of all elements in the system.
     *
     * @return (std::vector<std::string>) Element names.
     *
     * @code
     * auto names = engine.element_names();
     * std::cout << "Element 0: " << names[0] << std::endl;
     * @endcode
     */
    auto element_names() const -> std::vector<std::string>
    {
        return m_element_names;
    }

    /**
     * @brief Returns the names of all species in the system.
     *
     * @return (std::vector<std::string>) Species names.
     *
     * @code
     * auto names = engine.species_names();
     * std::cout << "Species 0: " << names[0] << std::endl;
     * @endcode
     */
    auto species_names() -> std::vector<std::string>
    {
        return m_species_names;
    }

    /**
     * @brief Returns the names of all phases in the system.
     *
     * @return (std::vector<std::string>) Phase names.
     *
     * @code
     * auto names = engine.phase_names();
     * std::cout << "Phase 0: " << names[0] << std::endl;
     * @endcode
     */
    auto phase_names() -> std::vector<std::string>
    {
        return m_phase_names;
    }

    /**
     * @brief Returns the aqueous phase name.
     *
     * @return (std::string) aqueous phase name. If empty, the aqueous phase is not in system.
     *
     * @code
     * auto aqueous_name = engine.aq_phase_symbol();
     * @endcode
    */
    auto aq_phase_symbol() -> std::string
    {
        return m_aq_phase_symbol;
    }

    /**
     * @brief Returns the gaseous phase name.
     *
     * @return (std::string) gaseous phase name. If empty, the gaseous phase is not in system.
     *
     * @code
     * auto gaseous_name = engine.gas_phase_symbol();
     * @endcode
    */
    auto gas_phase_symbol() -> std::string
    {
        return m_gas_phase_symbol;
    }

    /**
     * @brief Returns molar masses of elements.
     *
     * @return (ValuesMap) Dictionary of molar masses (kg/mol) for each element.
     *
     * @code
     * auto values = engine.element_molar_masses();
     * std::cout << "Molar mass of 'Fe': " << values["Fe"] << std::endl;
     * @endcode
     */
    auto element_molar_masses() const -> ValuesMap
    {
        return m_element_molar_masses;
    }

    /**
    * @brief Returns the names of all species for each phase in the system.
    *
    * @return (std::map<std::string, std::vector<std::string>>) Dictionary of species names.
    *
    * @code
    * auto map_names = engine.species_in_phase();
    * @endcode
    */
    auto species_in_phase() -> std::map<std::string, std::vector<std::string>>
    {
        return m_species_in_phase;
    }

    /**
     * @brief Returns the electrical charge of a species.
     *
    * @return (ValuesMap) Dictionary charges of the species.
     *
     * @code
     * auto charge = engine.species_charges();
     * std::cout << "Charge of 'OH-': " << charge["OH-"] << std::endl;
     * @endcode
     */
    auto species_charges() -> ValuesMap
    {
        return m_species_charges;
    }

    /**
     * @brief Returns molar masses of species.
     *
     * @return (ValuesMap) Dictionary of species molar masses (kg/mol).
     *
     * @code
    *  auto values = engine.species_molar_mass();
    *  std::cout << "Molar mass of 'H2O@': " << values["H2O@"] << std::endl;
     * @endcode
     */
    auto species_molar_mass() -> ValuesMap
    {
        return m_species_molar_mass;
    }

    /**
     * @brief Returns the standard molar volumes of a species.
     *
     * @return (ValuesMap) Dictionary of species standard molar volumes in m³/mol.
     *
     * @code
     * auto values = engine.species_molar_volumes();
     * std::cout << "Molar volume of 'H2O@': " << values["H2O@"] << std::endl;
     * @endcode
     */
    auto species_molar_volumes() -> ValuesMap
    {
        return m_species_molar_volumes;
    }

    /**
     * @brief Returns the amounts of the elements.
     *
     * @return (ValuesMap) Dictionary of elements amounts in mol.
     *
     * @code
     * auto amounts = engine.bulk_composition();
     * @endcode
     */
    auto bulk_composition() ->ValuesMap;

    /**
     * @brief Returns the pH of the aqueous phase.
     *
     * @return (double) pH (in the activity scale (-log10 molal)).
     *
     * @code
     * double ph = engine.pH();
     * @endcode
     */
    auto pH() -> double;

    /**
     * @brief Returns the pe of the aqueous phase.
     *
     * @return (double) pe (in the activity scale (-log10 molal)).
     *
     * @code
     * double pe = engine.pE();
     * @endcode
     */
    auto pE() -> double;

    /**
     * @brief Returns the ionic strength of the aqueous phase.
     *
     * @return (double) Ionic strength in molal.
     *
     * @code
     * double ionicStr = engine.ionic_strength();
     * @endcode
     */
    auto ionic_strength() -> double;

    /**
     * @brief Returns the total volume of the system.
     *
     * @return (double) Volume in m³.
     *
     * @code
     * double volume = engine.system_volume();
     * @endcode
     */
    auto system_volume() -> double;

    /**
     * @brief Returns the total mass of the system.
     *
     * @return (double) System mass in kg.
     *
     * @code
     * double mass = engine.system_mass();
     * @endcode
     */
    auto system_mass() -> double;

    /**
     * @brief Returns the Eh (redox potential) of the aqueous phase.
     *
     * @return (double) Eh in Volts.
     *
     * @code
     * double eh = engine.Eh();
     * @endcode
     */
    auto Eh() -> double;

    /**
     * @brief Checks if the last equilibrium computation converged.
     *
     * @return (bool) True if converged, false otherwise.
     *
     * @code
     * bool ok = engine.converged();
     * @endcode
     */
    auto converged() -> bool;

    /**
     * @brief Returns the number of iterations of the last equilibrium computation.
     *
     * @return (int) Number of iterations.
     *
     * @code
     * int n = engine.num_iterations();
     * @endcode
     */
    auto num_iterations() -> int;

    /**
     * @brief Returns the elapsed time of the last equilibrium computation.
     *
     * @return (double) Elapsed time in seconds.
     *
     * @code
     * double t = engine.elapsed_time();
     * @endcode
     */
    auto elapsed_time() -> double;

    /**
     * @brief Returns the total Gibbs energy of the system.
     *
     * @return (double) System Gibbs energy in J.
     *
     * @code
     * double G = engine.system_gibbs_energy();
     * @endcode
     */
    auto system_gibbs_energy() -> double;

    /**
     * @brief Returns the total enthalpy of the system.
     *
     * @return (double) System enthalpy in J.
     *
     * @code
     * double H = engine.system_enthalpy();
     * @endcode
     */
    auto system_enthalpy() -> double;

    /**
     * @brief Returns the total entropy of the system.
     *
     * @return (double) System entropy in J/K.
     *
     * @code
     * double S = engine.system_entropy();
     * @endcode
     */
    auto system_entropy() -> double;

    /**
     * @brief Returns the total isobaric heat capacity of the system.
     *
     * @return (double) System heat capacity in J/K.
     *
     * @code
     * double Cp = engine.system_heat_capacity_const_p();
     * @endcode
     */
    auto system_heat_capacity_const_p() -> double;

    /**
     * @brief Returns the per-element mass balance residual (mol) after equilibration.
     *
     * Computes W·n − b for each element, where W is the formula matrix, n is the
     * species amounts vector, and b is the stored input bulk composition. Values
     * near zero confirm a well-converged result; large values signal solver failure.
     *
     * @return (ValuesMap) Dictionary of element residuals in mol.
     *
     * @code
     * auto errors = engine.mass_balance_errors();
     * std::cout << "Residual for Ca: " << errors["Ca"] << " mol" << std::endl;
     * @endcode
     */
    auto mass_balance_errors() -> ValuesMap;

    /**
     * @brief Returns the per-element relative mass balance residual after equilibration.
     *
     * Computes (W·n − b) / b for each element. Elements whose input amount b ≈ 0
     * (e.g., the charge row Zz) are returned as 0.
     *
     * @return (ValuesMap) Dictionary of dimensionless relative residuals.
     *
     * @code
     * auto rel = engine.mass_balance_relative_errors();
     * std::cout << "Relative residual for Ca: " << rel["Ca"] << std::endl;
     * @endcode
     */
    auto mass_balance_relative_errors() -> ValuesMap;

    // -----------------------------------------------------------------------
    // Aqueous solution derived properties
    // -----------------------------------------------------------------------

    /**
     * @brief Returns the activity of water (H2O@) in the aqueous phase.
     *
     * @return (double) Water activity (dimensionless). Returns 1.0 if no aqueous phase.
     *
     * @code
     * double aw = engine.water_activity();
     * @endcode
     */
    auto water_activity() -> double;

    /**
     * @brief Returns the osmotic coefficient of the aqueous phase.
     *
     * Computed as Φ = −ln(a_w) × n_water / n_solutes.
     * Returns 1.0 for pure water or if no aqueous phase is present.
     *
     * @return (double) Osmotic coefficient (dimensionless).
     *
     * @code
     * double phi = engine.osmotic_coefficient();
     * @endcode
     */
    auto osmotic_coefficient() -> double;

    /**
     * @brief Returns the total hardness and Ca/Mg hardness components (mg/L as CaCO₃).
     *
     * Uses total Ca and Mg amounts in the aqueous phase (all species).
     * Returns keys: "total", "Ca", "Mg".
     *
     * @return (ValuesMap) Hardness dictionary in mg/L as CaCO₃.
     *
     * @code
     * auto hard = engine.hardness();
     * std::cout << "Total hardness: " << hard["total"] << " mg/L as CaCO3" << std::endl;
     * @endcode
     */
    auto hardness() -> ValuesMap;

    // -----------------------------------------------------------------------
    // Gas phase properties
    // -----------------------------------------------------------------------

    /**
     * @brief Returns partial pressures of all gas species (Pa).
     *
     * Computed as x_i · P_total for each species in the gas phase.
     * Returns an empty map if no gas phase is present.
     *
     * @return (ValuesMap) Dictionary of gas species partial pressures in Pa.
     *
     * @code
     * auto pp = engine.gas_species_partial_pressures();
     * std::cout << "CO2 partial pressure: " << pp["CO2"] << " Pa" << std::endl;
     * @endcode
     */
    auto gas_species_partial_pressures() -> ValuesMap;

    /**
     * @brief Returns fugacities of all gas species (Pa).
     *
     * Computed as exp(ln a_i) · P° where P° = 100000 Pa (1 bar standard state).
     * Returns an empty map if no gas phase is present.
     *
     * @return (ValuesMap) Dictionary of gas species fugacities in Pa.
     *
     * @code
     * auto f = engine.gas_species_fugacities();
     * std::cout << "CO2 fugacity: " << f["CO2"] << " Pa" << std::endl;
     * @endcode
     */
    auto gas_species_fugacities() -> ValuesMap;

    /**
     * @brief Returns fugacity coefficients of all gas species (dimensionless).
     *
     * Computed as exp(ln γ_i) where γ_i is the mole-fraction activity coefficient.
     * Returns an empty map if no gas phase is present.
     *
     * @return (ValuesMap) Dictionary of gas species fugacity coefficients.
     *
     * @code
     * auto phi = engine.gas_species_fugacity_coefficients();
     * std::cout << "CO2 fugacity coeff: " << phi["CO2"] << std::endl;
     * @endcode
     */
    auto gas_species_fugacity_coefficients() -> ValuesMap;

    /**
     * @brief Returns the partial pressure of a single gas species (Pa).
     *
     * @param species (std::string) Species name in the gas phase.
     * @return (double) Partial pressure in Pa, or 0 if species not in gas phase.
     *
     * @code
     * double pCO2 = engine.gas_partial_pressure("CO2");
     * @endcode
     */
    auto gas_partial_pressure(const std::string& species) -> double;

    /**
     * @brief Returns the fugacity of a single gas species (Pa).
     *
     * @param species (std::string) Species name in the gas phase.
     * @return (double) Fugacity in Pa, or 0 if species not in gas phase.
     *
     * @code
     * double fCO2 = engine.gas_fugacity("CO2");
     * @endcode
     */
    auto gas_fugacity(const std::string& species) -> double;

    /**
     * @brief Returns the fugacity coefficient of a single gas species (dimensionless).
     *
     * @param species (std::string) Species name in the gas phase.
     * @return (double) Fugacity coefficient, or 0 if species not in gas phase.
     *
     * @code
     * double phiCO2 = engine.gas_fugacity_coefficient("CO2");
     * @endcode
     */
    auto gas_fugacity_coefficient(const std::string& species) -> double;

    // -----------------------------------------------------------------------
    // Reaction thermodynamics
    // -----------------------------------------------------------------------

    /**
     * @brief Returns the standard log₁₀ equilibrium constant at current T and P.
     *
     * Computes log₁₀K = −ΔG°/(RT ln10) where ΔG° = Σ ν_i · G°_i.
     * The reaction dict maps species names to stoichiometric coefficients
     * (positive = product, negative = reactant).
     *
     * @param reaction (ValuesMap) Dictionary {species_name: stoich_coeff}.
     * @return (double) log₁₀K.
     *
     * @code
     * double logK = engine.log_K({{"Ca+2", 1}, {"CO3-2", 1}, {"Calcite", -1}});
     * @endcode
     */
    auto log_K(const ValuesMap& reaction) -> double;

    /**
     * @brief Returns the standard molar Gibbs energy change of a reaction at current T and P (J/mol).
     *
     * Computes ΔG° = Σ ν_i · G°_i.
     *
     * @param reaction (ValuesMap) Dictionary {species_name: stoich_coeff}.
     * @return (double) ΔG° in J/mol.
     *
     * @code
     * double dG = engine.delta_G0_reaction({{"Ca+2", 1}, {"CO3-2", 1}, {"Calcite", -1}});
     * @endcode
     */
    auto delta_G0_reaction(const ValuesMap& reaction) -> double;

    /**
     * @brief Returns the standard molar enthalpy change of a reaction at current T and P (J/mol).
     *
     * Computes ΔH° = Σ ν_i · H°_i.
     *
     * @param reaction (ValuesMap) Dictionary {species_name: stoich_coeff}.
     * @return (double) ΔH° in J/mol.
     *
     * @code
     * double dH = engine.delta_H0_reaction({{"Ca+2", 1}, {"CO3-2", 1}, {"Calcite", -1}});
     * @endcode
     */
    auto delta_H0_reaction(const ValuesMap& reaction) -> double;

    // -----------------------------------------------------------------------
    // Phase existence and type queries
    // -----------------------------------------------------------------------

    /**
     * @brief Returns true if a phase has a molar amount above the threshold.
     *
     * @param phase_name (std::string) Phase name.
     * @param threshold (double) Minimum amount in mol, default 1e-12.
     * @return (bool) True if the phase is present.
     *
     * @code
     * if (engine.phase_is_present("Calcite"))
     *     std::cout << "Calcite is precipitating" << std::endl;
     * @endcode
     */
    auto phase_is_present(const std::string& phase_name, double threshold = 1e-12) -> bool;

    /**
     * @brief Returns names of all phases whose molar amount exceeds the threshold.
     *
     * @param threshold (double) Minimum amount in mol, default 1e-12.
     * @return (std::vector<std::string>) Names of present phases.
     *
     * @code
     * for (auto& name : engine.present_phases())
     *     std::cout << name << std::endl;
     * @endcode
     */
    auto present_phases(double threshold = 1e-12) -> std::vector<std::string>;

    /**
     * @brief Returns names of all mineral (non-aqueous, non-gas) phases above the threshold.
     *
     * @param threshold (double) Minimum amount in mol, default 1e-12.
     * @return (std::vector<std::string>) Names of present mineral phases.
     *
     * @code
     * auto minerals = engine.present_minerals();
     * @endcode
     */
    auto present_minerals(double threshold = 1e-12) -> std::vector<std::string>;

    /**
     * @brief Returns true if the phase is the aqueous phase.
     *
     * @param phase (std::string) Phase name.
     * @return (bool) True if this is the aqueous phase.
     *
     * @code
     * engine.is_aqueous_phase("aq_gen"); // true
     * @endcode
     */
    auto is_aqueous_phase(const std::string& phase) -> bool
    {
        return !m_aq_phase_symbol.empty() && phase == m_aq_phase_symbol;
    }

    /**
     * @brief Returns true if the phase is the gas phase.
     *
     * @param phase (std::string) Phase name.
     * @return (bool) True if this is the gas phase.
     *
     * @code
     * engine.is_gas_phase("gas_gen"); // true
     * @endcode
     */
    auto is_gas_phase(const std::string& phase) -> bool
    {
        return !m_gas_phase_symbol.empty() && phase == m_gas_phase_symbol;
    }

    /**
     * @brief Returns true if the phase is a mineral (not aqueous, not gas).
     *
     * @param phase (std::string) Phase name.
     * @return (bool) True if this is a mineral phase.
     *
     * @code
     * engine.is_mineral_phase("Calcite"); // true
     * @endcode
     */
    auto is_mineral_phase(const std::string& phase) -> bool
    {
        return !is_aqueous_phase(phase) && !is_gas_phase(phase);
    }

    /**
     * @brief Returns the molar volumes of the phases.
     *
     * @return (ValuesMap) Dictionary of phases molar volumes in m³/mol.
     *
     * @code
     * auto volumes = engine.phases_molar_volume();
     * std::cout << "Molar volume of 'aq_gen': " << volumes["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_molar_volume() -> ValuesMap;

    /**
     * @brief Returns the saturation indices of all phases (log₁₀ units).
     *
     * @return (ValuesMap) Dictionary of phases saturation indices.
     *
     * @code
     * auto sat_indices = engine.phase_sat_indices();
     * std::cout << "Saturation indice of 'aq_gen': " << sat_indices["aq_gen"] << std::endl;
     * @endcode
     */
    auto phase_sat_indices() -> ValuesMap;

    /**
     * @brief Returns the densities of all phases (kg/m³).
     *
     * @return (ValuesMap) Dictionary of phase densities.
     *
     * @code
     * auto dens = engine.phases_density();
     * std::cout << "Density of 'aq_gen': " << dens["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_density() -> ValuesMap;

    /**
     * @brief Returns the total enthalpies of all phases (J).
     *
     * @return (ValuesMap) Dictionary of phase total enthalpies.
     *
     * @code
     * auto enth = engine.phases_enthalpy();
     * std::cout << "Enthalpy of 'aq_gen': " << enth["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_enthalpy() -> ValuesMap;

    /**
     * @brief Returns the total entropies of all phases (J/K).
     *
     * @return (ValuesMap) Dictionary of phase total entropies.
     *
     * @code
     * auto entr = engine.phases_entropy();
     * std::cout << "Entropy of 'aq_gen': " << entr["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_entropy() -> ValuesMap;

    /**
     * @brief Returns the molar Gibbs energies of all phases (J/mol).
     *
     * @return (ValuesMap) Dictionary of phase molar Gibbs energies.
     *
     * @code
     * auto gibbs = engine.phases_molar_gibbs_energy();
     * std::cout << "Molar Gibbs energy of 'aq_gen': " << gibbs["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_molar_gibbs_energy() -> ValuesMap;

    /**
     * @brief Returns the total isobaric heat capacities of all phases (J/K).
     *
     * @return (ValuesMap) Dictionary of phase total heat capacities.
     *
     * @code
     * auto cp = engine.phases_heat_capacity_const_p();
     * std::cout << "Heat capacity of 'aq_gen': " << cp["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_heat_capacity_const_p() -> ValuesMap;

     /**
     * @brief Returns the aq solution composition in mol/L aq solution.
     *
     * @return (ValuesMap) Dictionary for aq elements.
     *
     * @code
     * auto molarity = engine.aq_elements_molarity();
     * std::cout << "Molarity of 'O': " << molarity["O"] << std::endl;
     * @endcode
     */
    auto aq_elements_molarity() -> ValuesMap;

     /**
     * @brief Returns the aq solution elemental composition in mol/kgH2O.
     *
     * @return (ValuesMap) Dictionary for aq elements.
     *
     * @code
     * auto molality = engine.aq_elements_molality();
     * std::cout << "Molality of 'O': " << molality["O"] << std::endl;
     * @endcode
     */
    auto aq_elements_molality() -> ValuesMap;

    /**
     * @brief Returns the aq solution composition in mol/L of aqueous species.
     *
     * @return (ValuesMap) Dictionary for aq species.
     *
     * @code
     * auto molarity = engine.aq_species_molarity();
     * std::cout << "Molarity of 'H2O@': " << molarity["H2O@"] << std::endl;
     * @endcode
     */
    auto aq_species_molarity() -> ValuesMap;

    /**
     * @brief Returns the aq solution composition in mol/kg H2O of aqueous species (speciation).
     *
     * @return (ValuesMap) Dictionary for aq species.
     *
     * @code
     * auto molality = engine.aq_species_molality();
     * std::cout << "Molality of 'Mg(CO3)@': " << molality["Mg(CO3)@"] << std::endl;
     * @endcode
     */
    auto aq_species_molality() -> ValuesMap;

    /**
    * @brief Returns the amounts of elements in a aqueous phase.
    *
    * @return (ValuesMap) Dictionary containing the amounts of each element in the aqueous phase (in mol).
    *
    * @code
    * auto amount = engine.aq_elements_moles();
    * std::cout << "Amount of 'O' in the aqueous phase: " << amount["O"] << std::endl;
    * @endcode
    */
    auto aq_elements_moles() -> ValuesMap;

    /**
    * @brief Returns the mole amounts of elements in all solids together.
    *
    * @param min_amount_phase (double) min amount of phase in mol, default 1e-12.
    * @param min_amount_element (double) min amount of element in mol, default 1e-14.
    * @return (ValuesMap) Dictionary containing mole amounts of elements in all solids together.
    *
    * @code
    * auto amount = engine.solids_elements_moles();
    * std::cout << "Amount of 'Sn' in the solids phases: " << amount["Sn"] << std::endl;
    * @endcode
    */
    auto solids_elements_moles(double min_amount_phase=1e-12, double min_amount_element=1e-14) -> ValuesMap;

    /**
    * @brief Returns a dictionary (table) containing amounts of elements in phases in moles.
    *
    * @return (PhaseValuesMap) A dictionary of dictionaries containing mole amounts of elements for each phase.
    *
    * @code
    * auto phase_el_moles = engine.phases_elements_moles();
    * std::cout << "Amount of 'Ca' in 'aq_gen' phase: " << phase_el_moles["aq_gen"]["Ca"] << std::endl;
    * @endcode
    */
    auto phases_elements_moles() -> PhaseValuesMap;

    /**
     * @brief Returns the molar amounts of all phases.
     *
     * @return (ValuesMap) Dictionary of phases amounts in mol.
     *
     * @code
     * auto amounts = engine.phases_moles();
     * std::cout << "Amount of 'aq_gen' phase: " << amounts["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_moles() -> ValuesMap;

    /**
     * @brief Returns the amounts of all species.
     *
     * @return (ValuesMap) Dictionary of species amounts in mol.
     *
     * @code
     * auto amounts = engine.species_moles();
     * std::cout << "Amount of 'H2O@': " << amounts2["H2O@"] << std::endl;
     * @endcode
     */
    auto species_moles() -> ValuesMap;

    /**
     * @brief Returns the ln activities of all species.
     *
     * @return (ValuesMap) Dictionary of species ln Activities.
     *
     * @code
     * auto lnActivities = engine.species_ln_activities();
     * std::cout << "ln Activity of 'H2O@': " << lnActivities["H2O@"] << std::endl;
     * @endcode
     */
    auto species_ln_activities() -> ValuesMap;

    /**
     * @brief Returns the ln activity coefficients of all species (mole fraction scale).
     *
     * @return (ValuesMap) Dictionary of species ln Activity coefficients.
     *
     * @code
     * auto lnActCoeff = engine.species_ln_activity_coefficients();
     * std::cout << "ln Activity coeff of 'H2O@': " << lnActCoeff["H2O@"] << std::endl;
     * @endcode
     */
    auto species_ln_activity_coefficients() -> ValuesMap;

    /**
     * @brief Returns the mole fractions of all species.
     *
     * @return (ValuesMap) Dictionary of species mole fractions.
     *
     * @code
     * auto xfrac = engine.species_mole_fractions();
     * std::cout << "Mole fraction of 'H2O@': " << xfrac["H2O@"] << std::endl;
     * @endcode
     */
    auto species_mole_fractions() -> ValuesMap;

    /**
     * @brief Returns the ln concentrations of all species.
     *
     * @return (ValuesMap) Dictionary of species ln concentrations.
     *
     * @code
     * auto lnc = engine.species_ln_concentrations();
     * std::cout << "ln concentration of 'H2O@': " << lnc["H2O@"] << std::endl;
     * @endcode
     */
    auto species_ln_concentrations() -> ValuesMap;

    /**
     * @brief Returns the chemical potentials of all species (J/mol).
     *
     * @return (ValuesMap) Dictionary of species chemical potentials.
     *
     * @code
     * auto mu = engine.species_chemical_potentials();
     * std::cout << "Chemical potential of 'H2O@': " << mu["H2O@"] << std::endl;
     * @endcode
     */
    auto species_chemical_potentials() -> ValuesMap;

    /**
     * @brief Returns the standard molar Gibbs energies of all species (J/mol).
     *
     * @return (ValuesMap) Dictionary of standard molar Gibbs energies.
     *
     * @code
     * auto g0 = engine.species_gibbs_energies();
     * std::cout << "Standard molar Gibbs energy of 'H2O@': " << g0["H2O@"] << std::endl;
     * @endcode
     */
    auto species_gibbs_energies() -> ValuesMap;

    /**
     * @brief Returns the standard molar enthalpies of all species (J/mol).
     *
     * @return (ValuesMap) Dictionary of standard molar enthalpies.
     *
     * @code
     * auto h0 = engine.species_enthalpies();
     * std::cout << "Standard molar enthalpy of 'H2O@': " << h0["H2O@"] << std::endl;
     * @endcode
     */
    auto species_enthalpies() -> ValuesMap;

    /**
     * @brief Returns the standard molar entropies of all species (J/mol/K).
     *
     * @return (ValuesMap) Dictionary of standard molar entropies.
     *
     * @code
     * auto s0 = engine.species_entropies();
     * std::cout << "Standard molar entropy of 'H2O@': " << s0["H2O@"] << std::endl;
     * @endcode
     */
    auto species_entropies() -> ValuesMap;

    /**
     * @brief Returns the standard molar isobaric heat capacities of all species (J/mol/K).
     *
     * @return (ValuesMap) Dictionary of standard molar heat capacities.
     *
     * @code
     * auto cp0 = engine.species_heat_capacities_const_p();
     * std::cout << "Standard molar heat capacity of 'H2O@': " << cp0["H2O@"] << std::endl;
     * @endcode
     */
    auto species_heat_capacities_const_p() -> ValuesMap;

    /**
     * @brief Returns the upper limits for all species.
     *
     * @return (ValuesMap) Dictionary of species upper limits in mol.
     *
     * @code
     * auto upperLimits = engine.species_upper_bounds();
     * std::cout << "Upper bounds of 'Mg(HCO3)+': " << upperLimits["Mg(HCO3)+"] << std::endl;
     * @endcode
     */
    auto species_upper_bounds() -> ValuesMap;

    /**
     * @brief Returns the lower limits for all species.
     *
     * @return (ValuesMap) Dictionary of species lower limits in mol.
     *
     * @code
     * auto lowerLimits = engine.species_lower_bounds();
     * std::cout << "Lower bounds of 'Mg(HCO3)+': " << lowerLimits["Mg(HCO3)+"] << std::endl;
     * @endcode
     */
    auto species_lower_bounds() -> ValuesMap;

    /**
     * @brief Returns the amounts of the phase species.
     *
     * @param phase_symbol (std::string) phase name.
     * @return (ValuesMap) Dictionary of phase species amounts in mol.
     *
     * @code
     * auto amounts = engine.phase_species_moles("aq_gen");
     * std::cout << "Amount of 'H2O@' in 'aq_gen' phase:" << amounts["H2O@"] << std::endl;
     * @endcode
     */
    auto phase_species_moles(std::string phase_symbol) -> ValuesMap;

    /**
     * @brief Returns the mass(phase)/mass(system) ratios for [solid] phases.
     *
     * @return (ValuesMap) Dictionary of solids phases mass fraction.
     *
     * @code
     * auto mas_frac = engine.solids_mass_frac();
     * std::cout << "Mass fraction of 'Tin' phase: " << mas_frac["Tin"] << std::endl;
     * @endcode
     */
    auto solids_mass_frac() -> ValuesMap;

    /**
     * @brief Returns the volume(phase)/volume(total) ratio for solid phases.
     *
     * @return (ValuesMap) Dictionary of solids phases volume fraction from total system volume.
     *
     * @code
     * auto vol_frac = engine.solids_volume_frac();
     * std::cout << "Volume fraction of 'Tin' phase from total system volume: " << vol_frac["Tin"] << std::endl;
     * @endcode
     */
    auto solids_volume_frac() -> ValuesMap;

    /**
     * @brief Returns the volume fraction of aqueous phase from total system volume.
     *
     * @return (double) Volume fraction of aqueous phase from total system volume.
     *
     * @code
     * auto volume = engine.aq_volume_frac();
     * @endcode
     */
    auto aq_volume_frac() -> double;

    /**
     * @brief Returns the volumes of all phases.
     *
     * @return (ValuesMap) Dictionary of phases volumes in m³.
     *
     * @code
     * auto volumes = engine.phases_volume();
     * std::cout << "Volume of 'Tin' phase: " << volumes["Tin"] << std::endl;
     * @endcode
     */
    auto phases_volume() -> ValuesMap;

    /**
     * @brief Returns the masses of all phases.
     *
     * @return (ValuesMap) Dictionary of phases masses in kg.
     *
     * @code
     * auto masses = engine.phases_mass();
     * std::cout << "Mass of 'Tin' phase: " << masses["Tin"] << std::endl;
     * @endcode
     */
    auto phases_mass() -> ValuesMap;

    /**
     * @brief Returns the volume fractions of all phases from total system volume.
     *
     * @return (ValuesMap) Dictionary of phases and their volume fractions from total system volume.
     *
     * @code
     * auto volumes = engine.phases_volume_frac();
     * std::cout << "Volume fraction from total system volume of 'aq_gen' phase: " << volumes["aq_gen"] << std::endl;
     * @endcode
     */
    auto phases_volume_frac() -> ValuesMap;


    /**
     * @brief Add multiple species amounts in the system useful for adding aqueous solution composition.
     *
     * @param input_dict (ValuesMap) Dictionary of species amount in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     *
     * @code
     * engine.add_multiple_species_amt({ {"HCl@",0.01}, {"H2@",2} }, "moles");
     * @endcode
     */
    auto add_multiple_species_amt(const ValuesMap &input_dict, const std::string& units = "moles") -> void;

    /**
     * @brief Add species amount in the system useful for adding aqueous solution composition.
     *
     * @param species (std::string) Species symbol.
     * @param val (double) Species amount in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * engine.add_species_amt("H2O@", 0.01, "kg");
     * @endcode
     */
    auto add_species_amt(const std::string &species, double val, const std::string &units = "moles", std::optional<std::string> phase= std::nullopt) -> void;

    /**
     * @brief Add element amount in the system.
     *
     * @param element_name (std::string) Element symbol.
     * @param val (double) Element amount in units.
     * @param units (std::string) Units of amount ("moles", "kg"), default "moles".
     *
     * @code
     * engine.add_element_amt("Al", 0.3, "moles");
     * @endcode
     */
    auto add_element_amt(const std::string &element_name, double val, const std::string &units = "moles") -> void;

    /**
     * @brief Add multiple elements amount in the system useful for adding aqueous solution composition.
     *
     * @param input_dict (ValuesMap) Dictionary of elements amount in units.
     * @param units (std::string) Units of amount ("moles", "kg"), default "moles".
     *
     * @code
     * engine.add_multiple_elements_amt({ {"Na",1.013077}, {"Si",1.013077} }, "moles");
     * @endcode
     */
    auto add_multiple_elements_amt(const ValuesMap &input_dict, const std::string &units = "moles") -> void;

    /**
     * @brief Add multiple elements using user defined formula.
     *
     * @param formula (ValuesMap) User defined formula.
     * @param val (double) Component amount in units.
     * @param units (std::string) Units of amount ("moles", "kg"), default "moles".
     *
     * @code
     * engine.add_amt_from_formula( { {"K",2}, {"O",1} }, 4.108*1e-3, "kg");
     * @endcode
     */
    auto add_amt_from_formula(const ValuesMap &formula, double val, const std::string &units = "moles") -> void;

    /**
     * @brief Returns a bulk vector b from user-defined formula (as dict. {"H":2,"O":1} )
     * and amount of the formula [object] in units of 'moles' or 'kg'.
     *
     * @param formula (ValuesMap) User defined formula.
     * @param val (double) Amount of the formula [object] in units, default 1.
     * @param units (std::string) Units of amount ("moles", "kg"), default "moles".
     * @param min_amount (double) min amount of element in mol, default 1e-15.
     * @return (VectorConstRef) Vector of element amounts in mol.
     *
     * @code
     * auto vect = engine.get_b_from_formula( {{"H",2},{"O",1}}, 0.1, "kg");
     * @endcode
     */
    auto get_b_from_formula(const ValuesMap &formula, double val = 1, const std::string &units = "moles", double min_amount=1e-15) -> Vector;

    /**
     * @brief Sets an lower bound for multiple species.
     *
     * @param input_dict (ValuesMap) Dictionary of species lower bound.
     * @param units (std::string) Units of lower bound ("moles", "kg", "m3"), default "moles".
     *
     * @code
     * engine.set_multiple_species_lower_bound( {{"Mg(CO3)@",30}, {"Mg(HCO3)+",40}, {"Mg+2",50}});
     * @endcode
     */
    auto set_multiple_species_lower_bound(const ValuesMap &input_dict, const std::string &units = "moles") -> void;

    /**
     * @brief Sets an upper bounds for multiple species.
     *
     * @param input_dict (ValuesMap) Dictionary of species upper bound.
     * @param units (std::string) Units of upper bound ("moles", "kg", "m3"), default "moles".
     *
     * @code
     * engine.set_multiple_species_upper_bound( {{"Mg(CO3)@",300}, {"Mg(HCO3)+",400}, {"Mg+2",500}}, "moles");
     * @endcode
     */
    auto set_multiple_species_upper_bound(const ValuesMap &input_dict, const std::string &units = "moles") -> void;

    /**
     * @brief Sets a lower bound for a species identified by name.
     *
     * @param name (std::string) Species name.
     * @param val (double) Lower limit in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * engine.set_species_lower_bound( "Ca(HCO3)+", 200, "moles");
     * @endcode
     */
    auto set_species_lower_bound(const std::string& name, double val, const std::string& units= "moles",
                                 std::optional<std::string> phase= std::nullopt) -> void;

    /**
     * @brief Sets an upper bound for a species identified by name.
     *
     * @param name (std::string) Species name.
     * @param val (double) Upper limit in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * engine.set_species_upper_bound("CaOH+", 500, "kg");
     * @endcode
     */
    auto set_species_upper_bound(const std::string& name, double val, const std::string& units= "moles",
                                 std::optional<std::string> phase= std::nullopt) -> void;

    /**
     * @brief Sets a lower bound (minimum amount allowed to form) for a species identified by its index (phase depended case).
     *
     * @param ispecies (Index) Species index.
     * @param val (double) Lower limit in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     *
     * @code
     * engine.set_species_lower_bound(8, 400, "moles");
     * @endcode
     */
    auto set_species_lower_bound(Index ispecies, double val, const std::string& units= "moles") -> void;

    /**
     * @brief Sets an upper bound (maximum amount allowed to form) for a species identified by its index
     *  (phase depended case).
     *
     * @param ispecies (Index) Species index.
     * @param val (double) Upper limit in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     *
     * @code
     * engine.set_species_upper_bound( 8, 900, "kg");
     * @endcode
     */
     auto set_species_upper_bound(Index ispecies, double val, const std::string& units= "moles") -> void;


     /**
     * @brief Suppresses a phase in GEM calculation.
     *
     * @param phase_name (std::string) Phase name.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.suppress_phase("gas_gen");
     * @endcode
     */
    auto suppress_phase(const std::string &phase_name, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Suppresses multiple phases in calculation as given in phase names list.
     *
     * @param phase_name_list (std::vector<std::string>) Phases name list.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.suppress_multiple_phases({"Dolomite-dis", "Tin"});
     * @endcode
     */
    auto suppress_multiple_phases(const std::vector<std::string> &phase_name_list, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Suppresses a specie in calculation.
     *
     * @param species_name (std::string) Species name.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * engine.suppress_species("Ca(CO3)@");
     * @endcode
     */
    auto suppress_species(const std::string &species_name, double min_amount=0, double max_amount=1e-15,
                          std::optional<std::string> phase= std::nullopt) -> void;

    /**
     * @brief Suppresses multiple species in in GEM calculation as given in species name list.
     *
     * @param species_list (std::vector<std::string>) Species name list.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.suppress_multiple_species({"ClO4-", "Cl-"});
     * @endcode
     */
    auto suppress_multiple_species(const std::vector<std::string> &species_list, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Activate suppressed phase in GEM calculation.
     *
     * @param phase_name (std::string) Phase name.
     *
     * @code
     * engine.activate_phase("gas_gen");
     * @endcode
     */
    auto activate_phase(const std::string &phase_name) -> void;

    /**
     * @brief Activate multiple suppressed phases given in list.
     *
     * @param phase_name_list (std::vector<std::string>) Phases name list.
     *
     * @code
     * engine.activate_multiple_phases({"Dolomite-dis", "Tin"});
     * @endcode
     */
    auto activate_multiple_phases(const std::vector<std::string> &phase_name_list) -> void;

    /**
     * @brief Activate multiple suppressed species given in the list.
     *
     * @param species_list (std::vector<std::string>) Species name list.
     *
     * @code
     * engine.activate_multiple_species({"Ca(HCO3)+", "CaOH+", "Mg(CO3)@", "Mg(HCO3)+", "Mg+2", "ClO4-", "Cl-"});
     * @endcode
     */
    auto activate_multiple_species(const std::vector<std::string> &species_list) -> void;

    /**
     * @brief Activate a suppressed species in phase.
     *
     * @param species_name (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * engine.activate_species("Ca(CO3)@");
     * @endcode
     */
    auto activate_species(const std::string &species_name, std::optional<std::string> phase= std::nullopt) -> void;

    /**
     * @brief Returns all species amounts in moles.
     *
     * @return (PhaseValuesMap) A dictionary of dictionaries containing species amounts in mol for each phase.
     *
     * @code
     * auto amounts = engine.phase_species_moles();
     * std::cout << "Amount of 'CaOH+' in 'aq_gen' phase: " << amounts["aq_gen"]["CaOH+"] << std::endl;
     * @endcode
     */
    auto phase_species_moles() -> PhaseValuesMap;

    /**
     * @brief Returns the ln activities of the species.
     *
     * @return (PhaseValuesMap) A dictionary of dictionaries containing species ln Activities for each phase.
     *
     * @code
     * auto lnActivities = engine.phase_species_ln_activities();
     * std::cout << "ln activities of 'CaOH+' in 'aq_gen' phase: " << lnActivities["aq_gen"]["CaOH+"] << std::endl;
     * @endcode
     */
    auto phase_species_ln_activities() -> PhaseValuesMap;

    /**
     * @brief Returns the ln activity coefficients of the species (mole fraction scale).
     *
     * @return (PhaseValuesMap) A dictionary of dictionaries containing species ln Activity coefficients for each phase.
     *
     * @code
     * auto lnActCoeff = engine.phase_species_ln_activity_coefficients();
     * std::cout << "ln activities coeff of 'CaOH+' in 'aq_gen' phase: " << lnActCoeff["aq_gen"]["CaOH+"] << std::endl;
     * @endcode
     */
    auto phase_species_ln_activity_coefficients() -> PhaseValuesMap;

    /**
     * @brief Returns the upper limits for all species.
     *
     * @return (PhaseValuesMap) A dictionary of dictionaries containing species upper limits in mol for each phase.
     *
     * @code
     * auto upperLimits = engine.phase_species_upper_bounds();
     * std::cout << "Upper limits of 'CaOH+' in 'aq_gen' phase: " << upperLimits["aq_gen"]["CaOH+"] << std::endl;
     * @endcode
     */
    auto phase_species_upper_bounds() -> PhaseValuesMap;

    /**
     * @brief Returns the lower limits for all species.
     *
     * @return (PhaseValuesMap) A dictionary of dictionaries containing species lower limits in mol for each phase.
     *
     * @code
     * auto lowerLimits = engine.phase_species_lower_bounds();
     * std::cout << "Lower limits of 'CaOH+' in 'aq_gen' phase: " << lowerLimits["aq_gen"]["CaOH+"] << std::endl;
     * @endcode
     */
    auto phase_species_lower_bounds() -> PhaseValuesMap;

    friend auto operator<<(std::ostream &out, const ChemicalEngineMaps &state) -> std::ostream &;

protected:

    Vector b_amounts; // moles
    std::string input_file;
    mutable ChemicalEngine gem;

    std::vector<std::string> m_element_names;
    std::vector<std::string> m_species_names;
    std::vector<std::string> m_phase_names;

    std::string m_aq_phase_symbol;
    std::string m_gas_phase_symbol;

    /// dictionary containing elements and their molar masses
    ValuesMap m_element_molar_masses;
    /// dictionary containing species in phase
    std::map<std::string, std::vector<std::string>> m_species_in_phase;
    ///
    ValuesMap m_species_charges;
    /// dictionary containing species and their molar masses
    ValuesMap m_species_molar_mass;
    ///  dictionary containing species and their molar volumes
    ValuesMap m_species_molar_volumes;

    auto to_map( const std::vector<std::string>& names,  Vector values ) -> ValuesMap;

    auto to_phase_species_map( Vector values ) -> PhaseValuesMap;

    auto clear_vector(Vector& bb, double min_amount) -> void;

    auto init_caches(bool reset_calc, bool cold_start) -> void;
};

}
