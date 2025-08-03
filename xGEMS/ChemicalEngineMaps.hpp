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
 * https://bitbucket.org/gems4/xgems/src/master/demos/
 *
 *
 * @author R.A.Patel, Dmitrii Kulik, G.D. Miron, S.Dmytriieva
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
     * @param name (std::string) Species name.
     * @param value (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * engine.set_species_G0("H2O", -237140); // Value in J/mol.
     * @endcode
     */
    auto set_species_G0(std::string symbol, double value) -> void;

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
    auto element_names() -> std::vector<std::string>
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
    auto element_molar_masses() -> ValuesMap
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
     * @return (ValuesMap) Dictionary of solids phases volume fraction.
     *
     * @code
     * auto vol_frac = engine.solids_volume_frac();
     * std::cout << "Volume fraction of 'Tin' phase: " << vol_frac["Tin"] << std::endl;
     * @endcode
     */
    auto solids_volume_frac() -> ValuesMap;

    /**
     * @brief Returns the volume fraction of aqueous phase in the system.
     *
     * @return (double) Volume fraction of aqueous phase.
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
     * @brief Returns the volume fractions of all phases in the system.
     *
     * @return (ValuesMap) Dictionary of phases and their volume fractions .
     *
     * @code
     * auto volumes = engine.phases_volume_frac();
     * std::cout << "Volume fraction of 'aq_gen' phase: " << volumes["aq_gen"] << std::endl;
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
     *
     * @code
     * engine.add_species_amt("H2O@", 0.01, "kg");
     * @endcode
     */
    auto add_species_amt(const std::string &species, double val, const std::string &units = "moles") -> void;

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
     *
     * @code
     * engine.set_species_lower_bound( "Ca(HCO3)+", 200, "moles");
     * @endcode
     */
    auto set_species_lower_bound(const std::string& name, double val, const std::string& units= "moles") -> void;

    /**
     * @brief Sets an upper bound for a species identified by name.
     *
     * @param name (std::string) Species name.
     * @param amount (double) Upper limit in units.
     * @param units (std::string) Units of amount ("moles", "kg", "m3"), default "moles".
     *
     * @code
     * engine.set_species_upper_bound("CaOH+", 500, "kg");
     * @endcode
     */
    auto set_species_upper_bound(const std::string& name, double val, const std::string& units= "moles") -> void;

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
     * @brief Supresses a phase in GEM calculation.
     *
     * @param phase_name (std::string) Phase name.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.supress_phase("gas_gen");
     * @endcode
     */
    auto supress_phase(const std::string &phase_name, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Supresses multiple phases in calculation as given in phase names list.
     *
     * @param phase_name_list (std::vector<std::string>) Phases name list.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.supress_multiple_phases({"Dolomite-dis", "Tin"});
     * @endcode
     */
    auto supress_multiple_phases(const std::vector<std::string> &phase_name_list, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Supresses a specie in calculation.
     *
     * @param species_name (std::string) Species name.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.supress_species("Ca(CO3)@");
     * @endcode
     */
    auto supress_species(const std::string &species_name, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Supresses multiple species in in GEM calculation as given in species name list.
     *
     * @param species_list (std::vector<std::string>) Species name list.
     * @param min_amount (double) Lower amount of specie in mol, default 0.
     * @param max_amount (double) Upper amount of specie in mol, default 1e-15.
     *
     * @code
     * engine.supress_multiple_species({"ClO4-", "Cl-"});
     * @endcode
     */
    auto supress_multiple_species(const std::vector<std::string> &species_list, double min_amount=0, double max_amount=1e-15) -> void;

    /**
     * @brief Activate supressed phase in GEM calculation.
     *
     * @param phase_name (std::string) Phase name.
     *
     * @code
     * engine.activate_phase("gas_gen");
     * @endcode
     */
    auto activate_phase(const std::string &phase_name) -> void;

    /**
     * @brief Activate multiple supressed phases given in list.
     *
     * @param phase_name_list (std::vector<std::string>) Phases name list.
     *
     * @code
     * engine.activate_multiple_phases({"Dolomite-dis", "Tin"});
     * @endcode
     */
    auto activate_multiple_phases(const std::vector<std::string> &phase_name_list) -> void;

    /**
     * @brief Activate multiple supressed species given in the list.
     *
     * @param species_list (std::vector<std::string>) Species name list.
     *
     * @code
     * engine.activate_multiple_species({"Ca(HCO3)+", "CaOH+", "Mg(CO3)@", "Mg(HCO3)+", "Mg+2", "ClO4-", "Cl-"});
     * @endcode
     */
    auto activate_multiple_species(const std::vector<std::string> &species_list) -> void;

    /**
     * @brief Activate a supressed species in phase.
     *
     * @param species_name (std::string) Species name.
     *
     * @code
     * engine.activate_species("Ca(CO3)@");
     * @endcode
     */
    auto activate_species(const std::string &species_name) -> void;

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

protected:

    Vector b_amounts; // moles
    std::string input_file;
    ChemicalEngine gem;

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
};

}
