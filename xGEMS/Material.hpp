/**
 * @file Material.hpp
 * @brief Header file for the Material class — a user-friendly recipe builder
 *        for assembling bulk-composition (b) vectors used by ChemicalEngine.
 *
 * Build a recipe by chaining .add() calls with chemical formulas (or
 * element dictionaries) and amounts. Scale the whole recipe by calling
 * the material with a target total: ``material(1.0, "kg")``, or
 * rescale in place with ``setQuantity(1.0, "kg")``. Combine two
 * recipes with ``operator+``.
 *
 * Units follow the xGEMS convention used elsewhere in the library
 * (kg and mol internally). Supported strings:
 *   - Mass : "kg", "g", "mg", "ug"
 *   - Moles: "mol" (or "moles"), "mmol", "umol"
 *
 * @code
 *   xGEMS::ChemicalEngine engine("system-dat.lst");
 *
 *   // ── Seawater (Millero et al.) — H2O separated from dissolved ions ──
 *   xGEMS::Material seawater(engine, "seawater");
 *   seawater.add("H2O", 1.00, "kg")
 *           .add({{"Cl",0.5467499},{"Na",0.4689999},
 *                 {"Mg",0.05282088},{"Ca",0.01037339},
 *                 {"K", 0.01020831},{"S", 0.02823655},
 *                 {"C", 0.001966361}},
 *                "mol");   // absolute mol per element — no scale factor
 *
 *   // ── Basalt glass — oxide blend ──
 *   xGEMS::Material basalt(engine, "basalt_glass");
 *   basalt.add("SiO2",  1.000, "mol")
 *         .add("Al2O3", 0.175, "mol")
 *         .add("FeO",   0.169, "mol")
 *         .add("CaO",   0.270, "mol")
 *         .add("MgO",   0.260, "mol")
 *         .add("Na2O",  0.040, "mol")
 *         .add("K2O",   0.004, "mol")
 *         .add("TiO2",  0.020, "mol")
 *         .add("Fe2O3", 0.006, "mol");
 *
 *   // ── Leucogranite — mineral assemblage by mass ──
 *   xGEMS::Material granite(engine, "leucogranite");
 *   granite.add("KAlSi3O8",        17.0, "g")    // Microcline
 *          .add("KAl3Si3O10(OH)2", 19.0, "g")    // Muscovite
 *          .add("NaAlSi3O8",       29.0, "g")    // Albite
 *          .add("SiO2",            35.0, "g");   // Quartz
 *
 *   // ── Combine and scale ──
 *   auto mix = seawater(1.0, "kg") + basalt(10.0, "g");
 *   engine.equilibrate(298.15, 1.0e5, mix.b());
 *
 *   // Rescale in place — useful when iterating on a target total.
 *   granite.setQuantity(0.5, "kg");
 * @endcode
 *
 * @author G.D. Miron and contributors
 * @date 2026
 *
 * license GNU General Public License v3 or later
 */

#pragma once

// C++ includes
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

// xGEMS includes
#include <xGEMS/ChemicalEngine.hpp>
#include <xGEMS/ChemicalEngineMaps.hpp>
#include <xGEMS/Eigen.hpp>

namespace xGEMS
{

/**
 * @brief Alias for {element_name: amount} maps used as the natural Python
 *        analogue for element dictionaries.
 */
using ElementMap = std::map<std::string, double>;


/**
 * @class MaterialEngineAdapter
 * @brief Polymorphic adapter that lets Material work with either
 *        ChemicalEngine or ChemicalEngineMaps uniformly.
 *
 * Concrete subclasses are constructed internally by Material and shared
 * across copies via ``std::shared_ptr``. Users typically don't see this
 * type — they just pass an engine to the Material constructor.
 *
 * Element names and molar masses are cached at construction so subsequent
 * lookups are O(log n) regardless of which engine variant is wrapped.
 */
class MaterialEngineAdapter
{
public:
    virtual ~MaterialEngineAdapter() = default;

    /// Element names in the engine's index order.
    virtual auto elementNames() const -> const std::vector<std::string>& = 0;

    /// Molar mass [kg/mol] for one element, looked up by name.
    /// Throws ``std::out_of_range`` if the name isn't in the engine.
    virtual auto elementMolarMass(const std::string& name) const -> double = 0;

    /// Index of an element by name, or -1 if absent.
    virtual auto indexOf(const std::string& name) const -> Index = 0;
};

/**
 * @brief Parse a chemical formula string into element counts.
 *
 * Supports nested parentheses and hydrate notation with ``.`` separators
 * (e.g. ``"CaSO4.2H2O"``). Element symbols follow the standard pattern
 * uppercase letter, optional lowercase letter, optional integer subscript.
 *
 * @param formula The formula to parse (e.g. ``"Ca(OH)2"``).
 * @return (ElementMap) A map of element symbol to count (mol per
 *         formula unit).
 *
 * @code
 *   auto el = xGEMS::parseChemicalFormula("CaSO4.2H2O");
 *   // el == { {"Ca",1}, {"S",1}, {"O",6}, {"H",4} }
 * @endcode
 */
auto parseChemicalFormula(const std::string& formula) -> ElementMap;


/**
 * @class Material
 * @brief A recipe of chemical constituents bound to a ChemicalEngine.
 *
 * The Material stores an ordered list of constituents, each described
 * by a chemical formula (or an explicit element map), a value, and a
 * unit string. ``elements()`` and ``b()`` collapse the recipe into an
 * element-mole map / Eigen vector ready to be passed to
 * ``ChemicalEngine::equilibrate``.
 *
 * Scaling and combination follow the natural intuition:
 *   - ``material(value, units)``  → returns a scaled COPY whose total
 *                                   equals ``value`` in ``units``.
 *   - ``material.scale(factor)``  → in-place multiplicative rescaling.
 *   - ``a + b``                   → returns a new Material containing
 *                                   the union of both recipes.
 *
 * The engine reference is stored as a raw pointer; the engine must
 * outlive the Material(s) bound to it. A default-constructed Material
 * has no engine — useful as a return-value sink, but ``.b()`` will
 * throw until the engine is set.
 */
class Material
{
public:
    /**
     * @brief Default-construct an empty Material with no engine bound.
     *
     * ``.b()`` and any methods that depend on element-name lookup will
     * throw until the Material is given an engine (typically by being
     * the result of ``operator+`` or copy-assignment from a bound one).
     */
    Material();

    /**
     * @brief Construct an empty Material bound to a ChemicalEngine.
     *
     * @param engine The ChemicalEngine providing the element list and
     *               molar masses. Must outlive the Material.
     * @param name   Optional name used in ``toString()`` and combined
     *               material labels.
     *
     * @code
     *   xGEMS::ChemicalEngine engine("system-dat.lst");
     *   xGEMS::Material brine(engine, "brine");
     * @endcode
     */
    explicit Material(const ChemicalEngine& engine,
                      std::string name = "material");

    /**
     * @brief Construct an empty Material bound to a ChemicalEngineMaps.
     *
     * Same semantics as the ChemicalEngine overload; provided so the
     * dict-style API works seamlessly with the same Material class.
     *
     * @code
     *   xGEMS::ChemicalEngineMaps engine("system-dat.lst");
     *   xGEMS::Material brine(engine, "brine");
     * @endcode
     */
    explicit Material(const ChemicalEngineMaps& engine,
                      std::string name = "material");

    /**
     * @brief Add a constituent by chemical formula string.
     *
     * @param formula   Chemical formula (e.g. ``"CaCO3"``, ``"H2O"``,
     *                  ``"Cl"``, ``"CaSO4.2H2O"``).
     * @param value     Amount in ``units``; default 1.0.
     * @param units     Unit string. Default ``"mol"``. See class header.
     * @param react_ext Reaction-extent multiplier in [0, 1]; default 1.
     * @return (Material&) ``*this`` for chaining.
     *
     * @code
     *   brine.add("H2O", 1.00, "kg")
     *        .add("NaCl", 1.00, "mol");
     * @endcode
     */
    auto add(const std::string& formula,
             double value = 1.0,
             const std::string& units = "mol",
             double react_ext = 1.0) -> Material&;

    /**
     * @brief Add a constituent by an explicit element-count map.
     *
     * Equivalent to adding ``"NaCl"`` if you pass ``{{"Na",1},{"Cl",1}}``.
     * Useful when you want to specify a "synthetic recipe" whose element
     * proportions don't follow a clean chemical formula (e.g. a sea-salt
     * blend with mixed cation ratios).
     *
     * @param formula   Map of element symbol to count.
     * @param value     Amount in ``units``; default 1.0.
     * @param units     Unit string. Default ``"mol"``.
     * @param react_ext Reaction-extent multiplier in [0, 1]; default 1.
     * @return (Material&) ``*this`` for chaining.
     *
     * @code
     *   brine.add( { {"Na",1}, {"Cl",1} }, 1.0, "mol");
     * @endcode
     */
    auto add(const ElementMap& formula,
             double value = 1.0,
             const std::string& units = "mol",
             double react_ext = 1.0) -> Material&;

    /**
     * @brief Add elements from an absolute-amount map.
     *
     * When no scale value is given the dict values ARE the absolute
     * amounts of each element in ``units``. Each entry is stored as its
     * own single-element constituent; the total is the sum of all values
     * converted through ``units``.
     *
     * This is the natural form for specifying dissolved-ion concentrations
     * when the amounts are already known per-element:
     *
     * @param amounts   Map of element symbol to absolute amount in ``units``.
     * @param units     Unit string. Default ``"mol"``.
     * @param react_ext Reaction-extent multiplier in [0, 1]; default 1.
     * @return (Material&) ``*this`` for chaining.
     *
     * @code
     *   // Millero seawater dissolved ions — values are mol per kg H2O
     *   seawater.add({{"Cl",0.5467499},{"Na",0.4689999},
     *                 {"Mg",0.05282088},{"Ca",0.01037339},
     *                 {"K", 0.01020831},{"S", 0.02823655},
     *                 {"C", 0.001966361}},
     *                "mol");
     * @endcode
     */
    auto add(const ElementMap& amounts,
             const std::string& units,
             double react_ext = 1.0) -> Material&;

    /**
     * @brief Return a scaled COPY whose total equals ``(value, units)``.
     *
     * Computes the current total in the requested units, then
     * proportionally rescales every constituent's value. Useful for
     * the natural-looking ``brine(1.0, "kg") + rock(10.0, "g")`` syntax.
     *
     * @param value Target total amount; default 1.0.
     * @param units Target unit; default ``"mol"``.
     * @return (Material) Scaled copy. The original is unchanged.
     *
     * @throws std::runtime_error if the current total is non-positive.
     */
    auto operator()(double value = 1.0,
                    const std::string& units = "mol") const -> Material;

    /**
     * @brief Concatenate two recipes into a new Material.
     *
     * The result owns a copy of both constituent lists; the engine
     * pointer is inherited from whichever operand has one (left-hand
     * side takes precedence).
     */
    auto operator+(const Material& other) const -> Material;

    /**
     * @brief Multiply every constituent's value by ``factor`` in place.
     * @return (Material&) ``*this`` for chaining.
     */
    auto scale(double factor) -> Material&;

    /**
     * @brief Rescale in place so the total mass equals ``kg`` kilograms.
     * @return (Material&) ``*this`` for chaining.
     * @throws std::runtime_error if current total mass is non-positive.
     */
    auto scaleToMass(double kg) -> Material&;

    /**
     * @brief Rescale in place so the total amount equals ``mol`` moles.
     * @return (Material&) ``*this`` for chaining.
     * @throws std::runtime_error if current total amount is non-positive.
     */
    auto scaleToMoles(double mol) -> Material&;

    /**
     * @brief Rescale in place so the total equals ``(value, units)``.
     *
     * Companion to ``operator()(value, units)`` (which returns a scaled
     * copy). ``setQuantity`` modifies the Material directly and returns
     * ``*this`` for chaining — useful when iterating on a target total
     * or rescaling within a larger workflow.
     *
     * @param value Target total amount.
     * @param units Target unit; default ``"mol"``. Same set as elsewhere.
     * @return (Material&) ``*this`` for chaining.
     * @throws std::runtime_error if the current total is non-positive.
     */
    auto setQuantity(double value, const std::string& units = "mol") -> Material&;

    /**
     * @brief Return the recipe's contribution to the bulk composition
     *        as an element-mole map.
     *
     * Each constituent is resolved to (formula moles) × (element
     * counts) and summed across the recipe. Independent of the engine,
     * but if no engine is bound the molar-mass conversions used by
     * mass-based constituents will already have happened at .add() time.
     *
     * @return (ElementMap) {element_name: total moles in this material}
     */
    auto elements() const -> ElementMap;

    /**
     * @brief Return the bulk composition vector ``b`` indexed by the
     *        engine's element list.
     *
     * The returned vector has length ``engine.numElements()`` and is
     * ready to be passed to ``ChemicalEngine::equilibrate(T, P, b)``.
     *
     * @return (Vector) Element-mole vector matching engine indexing.
     * @throws std::runtime_error if no engine is bound, or if the
     *         recipe contains an element not in the engine's list.
     */
    auto b() const -> Vector;

    /**
     * @brief Return the bulk composition as a {element_name: mol} map.
     *
     * Natural form for ChemicalEngineMaps users. Includes one entry
     * per engine element (zeros for elements not in the recipe).
     *
     * @return (ElementMap) {element_name: total moles}.
     * @throws std::runtime_error if no engine is bound.
     */
    auto bMap() const -> ElementMap;

    /**
     * @brief Total mass of the recipe in kg.
     */
    auto massKg() const -> double;

    /**
     * @brief Total amount of the recipe in mol (sum of constituent moles).
     */
    auto moles() const -> double;

    /**
     * @brief Total in any supported unit.
     *
     * @param units Unit string (default ``"mol"``).
     */
    auto total(const std::string& units = "mol") const -> double;

    /**
     * @brief Number of constituents in the recipe.
     */
    auto size() const -> std::size_t;

    /**
     * @brief Default amount [mol] used as the initial value for every
     *        engine element in ``b()`` and ``bMap()``.
     *
     * GEMS requires all entries of the bulk-composition vector to be
     * strictly positive. Setting a small non-zero floor (default 1e-15)
     * ensures elements absent from the recipe still satisfy this constraint
     * without meaningfully affecting the equilibrium result.
     */
    auto defaultAmount() const -> double;
    auto setDefaultAmount(double amount) -> Material&;

    /**
     * @brief Material name (used in ``toString()`` and combined labels).
     */
    auto name() const -> const std::string&;
    auto setName(std::string n) -> Material&;

    /**
     * @brief Pretty one-line summary, e.g.
     *        ``"Material('brine', 5 constituents, 1078.7 g, 56.76 mol)"``.
     */
    auto toString() const -> std::string;

private:
    /**
     * @brief One entry in the recipe (formula + value + unit + react_ext).
     */
    struct Constituent
    {
        std::string name;          ///< Display name (formula string or "elements").
        ElementMap  formula;       ///< Element counts per 1 mol of formula.
        double      formula_mass_kg; ///< Molar mass (kg/mol), computed once.
        double      value;         ///< User-supplied value.
        std::string units;         ///< User-supplied unit string.
        double      react_ext;     ///< Reaction-extent multiplier in [0, 1].

        auto moles() const -> double;
        auto massKg() const -> double;
        auto elementMoles() const -> ElementMap;
    };

    const MaterialEngineAdapter* adapter() const { return adapter_.get(); }
    std::shared_ptr<const MaterialEngineAdapter> adapter_;
    std::string           name_{"material"};
    std::vector<Constituent> constituents_;
    double                default_amount_{1.0e-15};

    /// Compute the molar mass of a formula (kg/mol) using the engine's
    /// element-molar-mass table. Throws if no engine bound or if the
    /// formula refers to elements absent from the engine.
    auto formulaMassKg(const ElementMap& formula) const -> double;

    /// Build a constituent and validate its unit string.
    auto buildConstituent(const ElementMap& formula,
                          const std::string& name,
                          double value,
                          const std::string& units,
                          double react_ext) const -> Constituent;
};

} // namespace xGEMS
