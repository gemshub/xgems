/**
 * @file Material.cpp
 * @brief Implementation of the Material class — a recipe builder for
 *        ChemicalEngine bulk-composition (b) vectors.
 *
 * Internal units are kg and mol, matching the rest of the xGEMS API.
 * The formula parser is a small recursive-descent over uppercase-lowercase
 * element symbols, optional integer subscripts, nested parentheses, and
 * ``.``-separated hydrate parts.
 *
 * license GNU General Public License v3 or later
 */

#include "Material.hpp"

#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace xGEMS
{

// ─── Anonymous-namespace helpers ─────────────────────────────────────

namespace
{

// Conversion to canonical units (kg for mass, mol for amount).
const std::unordered_map<std::string, double> MASS_TO_KG = {
    {"kg", 1.0e+0 },
    {"g",  1.0e-3 },
    {"mg", 1.0e-6 },
    {"ug", 1.0e-9 },
    {"mkg",1.0e-9 },   // micro-kg ≡ μg, matches the GEMS UI dropdown
};

const std::unordered_map<std::string, double> MOLE_TO_MOL = {
    {"mol",   1.0e+0 },
    {"moles", 1.0e+0 },   // accepted synonym (matches ChemicalEngineMaps API)
    {"mmol",  1.0e-3 },
    {"umol",  1.0e-6 },
};

inline auto isMassUnit(const std::string& u) -> bool
{
    return MASS_TO_KG.find(u) != MASS_TO_KG.end();
}

inline auto isMoleUnit(const std::string& u) -> bool
{
    return MOLE_TO_MOL.find(u) != MOLE_TO_MOL.end();
}

inline auto checkUnit(const std::string& u) -> void
{
    if (!isMassUnit(u) && !isMoleUnit(u)) {
        throw std::invalid_argument(
            "xGEMS::Material: unsupported unit '" + u + "'. "
            "Use one of: mol, moles, mmol, umol, kg, g, mg, ug.");
    }
}

// ─── Concrete engine adapters ────────────────────────────────────────

/// Adapter wrapping the index-based ChemicalEngine. Caches the element
/// name list and a {name → kg/mol, index} side-table at construction.
class ChemEngineAdapter final : public MaterialEngineAdapter
{
public:
    explicit ChemEngineAdapter(const ChemicalEngine& engine)
    {
        const Index n = engine.numElements();
        names_.reserve(static_cast<std::size_t>(n));
        const auto mm = engine.elementMolarMasses();
        for (Index i = 0; i < n; ++i) {
            std::string nm = engine.elementName(i);
            names_.push_back(nm);
            molar_mass_[nm] = mm[i];
            index_[nm] = i;
        }
    }

    auto elementNames() const
        -> const std::vector<std::string>& override { return names_; }

    auto elementMolarMass(const std::string& name) const -> double override
    {
        auto it = molar_mass_.find(name);
        if (it == molar_mass_.end()) {
            throw std::runtime_error(
                "xGEMS::Material: element '" + name +
                "' not found in ChemicalEngine.");
        }
        return it->second;
    }

    auto indexOf(const std::string& name) const -> Index override
    {
        auto it = index_.find(name);
        return (it == index_.end()) ? Index{-1} : it->second;
    }

private:
    std::vector<std::string>      names_;
    std::map<std::string, double> molar_mass_;
    std::map<std::string, Index>  index_;
};


/// Adapter wrapping the dict-based ChemicalEngineMaps. Same cached
/// representation as the ChemicalEngine adapter — translation happens
/// once at construction.
class ChemEngineMapsAdapter final : public MaterialEngineAdapter
{
public:
    explicit ChemEngineMapsAdapter(const ChemicalEngineMaps& engine)
    {
        names_ = engine.element_names();
        auto mm = engine.element_molar_masses();   // ValuesMap (kg/mol)
        for (std::size_t i = 0; i < names_.size(); ++i) {
            const std::string& nm = names_[i];
            molar_mass_[nm] = mm.at(nm);
            index_[nm] = static_cast<Index>(i);
        }
    }

    auto elementNames() const
        -> const std::vector<std::string>& override { return names_; }

    auto elementMolarMass(const std::string& name) const -> double override
    {
        auto it = molar_mass_.find(name);
        if (it == molar_mass_.end()) {
            throw std::runtime_error(
                "xGEMS::Material: element '" + name +
                "' not found in ChemicalEngineMaps.");
        }
        return it->second;
    }

    auto indexOf(const std::string& name) const -> Index override
    {
        auto it = index_.find(name);
        return (it == index_.end()) ? Index{-1} : it->second;
    }

private:
    std::vector<std::string>      names_;
    std::map<std::string, double> molar_mass_;
    std::map<std::string, Index>  index_;
};

/// Recursive-descent parser for an atomic formula segment (no '.' separators).
/// Handles uppercase-lowercase element symbols, integer subscripts, and
/// nested parentheses.
auto parseSimple(const std::string& s) -> ElementMap
{
    std::vector<ElementMap> stack;
    stack.emplace_back();

    const std::size_t n = s.size();
    std::size_t i = 0;
    while (i < n) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == '(') {
            stack.emplace_back();
            ++i;
        }
        else if (c == ')') {
            if (stack.size() == 1)
                throw std::invalid_argument(
                    "xGEMS::Material: unmatched ')' in formula '" + s + "'.");
            ElementMap top = std::move(stack.back());
            stack.pop_back();
            ++i;
            int mult = 0;
            while (i < n && std::isdigit(static_cast<unsigned char>(s[i]))) {
                mult = mult * 10 + (s[i] - '0');
                ++i;
            }
            if (mult == 0) mult = 1;
            for (const auto& kv : top) {
                stack.back()[kv.first] += kv.second * mult;
            }
        }
        else if (std::isupper(c)) {
            std::string el(1, static_cast<char>(c));
            ++i;
            if (i < n && std::islower(static_cast<unsigned char>(s[i]))) {
                el += s[i];
                ++i;
            }
            int cnt = 0;
            while (i < n && std::isdigit(static_cast<unsigned char>(s[i]))) {
                cnt = cnt * 10 + (s[i] - '0');
                ++i;
            }
            if (cnt == 0) cnt = 1;
            stack.back()[el] += cnt;
        }
        else {
            throw std::invalid_argument(
                "xGEMS::Material: unsupported token in formula '" + s + "'.");
        }
    }
    if (stack.size() != 1)
        throw std::invalid_argument(
            "xGEMS::Material: unmatched '(' in formula '" + s + "'.");
    return stack.back();
}

} // anonymous namespace

// ─── parseChemicalFormula ────────────────────────────────────────────

auto parseChemicalFormula(const std::string& formula) -> ElementMap
{
    // Strip whitespace.
    std::string s;
    s.reserve(formula.size());
    for (char c : formula) {
        if (!std::isspace(static_cast<unsigned char>(c))) s.push_back(c);
    }

    // Split on '.' for hydrate notation.
    std::vector<std::string> parts;
    {
        std::string cur;
        for (char c : s) {
            if (c == '.') {
                parts.push_back(std::move(cur));
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
        parts.push_back(std::move(cur));
    }

    ElementMap total;
    for (const std::string& part : parts) {
        if (part.empty()) continue;
        // Optional leading integer coefficient ("2H2O" → coef = 2).
        std::size_t k = 0;
        while (k < part.size() &&
               std::isdigit(static_cast<unsigned char>(part[k]))) ++k;
        const double coef = (k > 0) ? std::stod(part.substr(0, k)) : 1.0;
        const ElementMap sub = parseSimple(part.substr(k));
        for (const auto& kv : sub) {
            total[kv.first] += kv.second * coef;
        }
    }
    return total;
}


// ─── Material::Constituent ───────────────────────────────────────────

auto Material::Constituent::moles() const -> double
{
    if (auto it = MOLE_TO_MOL.find(units); it != MOLE_TO_MOL.end()) {
        return value * it->second * react_ext;
    }
    if (auto it = MASS_TO_KG.find(units); it != MASS_TO_KG.end()) {
        return value * it->second / formula_mass_kg * react_ext;
    }
    throw std::logic_error(
        "xGEMS::Material::Constituent: invalid unit '" + units + "'.");
}

auto Material::Constituent::massKg() const -> double
{
    return moles() * formula_mass_kg;
}

auto Material::Constituent::elementMoles() const -> ElementMap
{
    const double mol = moles();
    ElementMap out;
    for (const auto& kv : formula) out[kv.first] = kv.second * mol;
    return out;
}


// ─── Material ────────────────────────────────────────────────────────

Material::Material() = default;

Material::Material(const ChemicalEngine& engine, std::string name)
  : adapter_(std::make_shared<ChemEngineAdapter>(engine)),
    name_(std::move(name))
{}

Material::Material(const ChemicalEngineMaps& engine, std::string name)
  : adapter_(std::make_shared<ChemEngineMapsAdapter>(engine)),
    name_(std::move(name))
{}

auto Material::name() const -> const std::string& { return name_; }
auto Material::setName(std::string n) -> Material& {
    name_ = std::move(n);
    return *this;
}
auto Material::size() const -> std::size_t { return constituents_.size(); }
auto Material::defaultAmount() const -> double { return default_amount_; }
auto Material::setDefaultAmount(double amount) -> Material& {
    default_amount_ = amount;
    return *this;
}


auto Material::formulaMassKg(const ElementMap& formula) const -> double
{
    if (!adapter_) {
        throw std::runtime_error(
            "xGEMS::Material: cannot compute molar mass without an engine. "
            "Construct with Material(engine).");
    }
    double m_kg = 0.0;
    for (const auto& kv : formula) {
        m_kg += kv.second * adapter_->elementMolarMass(kv.first);
    }
    return m_kg;
}

auto Material::buildConstituent(const ElementMap& formula,
                                const std::string& nm,
                                double value,
                                const std::string& units,
                                double react_ext) const -> Constituent
{
    checkUnit(units);
    Constituent c;
    c.name            = nm;
    c.formula         = formula;
    c.formula_mass_kg = formulaMassKg(formula);
    if (formula.empty() || c.formula_mass_kg <= 0.0) {
        throw std::invalid_argument(
            "xGEMS::Material: formula must contain at least one valid element.");
    }
    c.value     = value;
    c.units     = units;
    c.react_ext = react_ext;
    return c;
}

auto Material::add(const std::string& formula, double value,
                   const std::string& units, double react_ext) -> Material&
{
    auto parsed = parseChemicalFormula(formula);
    constituents_.push_back(
        buildConstituent(parsed, formula, value, units, react_ext));
    return *this;
}

auto Material::add(const ElementMap& formula, double value,
                   const std::string& units, double react_ext) -> Material&
{
    constituents_.push_back(
        buildConstituent(formula, "elements", value, units, react_ext));
    return *this;
}

auto Material::add(const ElementMap& amounts, const std::string& units,
                   double react_ext) -> Material&
{
    checkUnit(units);
    for (const auto& [el, amt] : amounts) {
        constituents_.push_back(
            buildConstituent({{el, 1.0}}, el, amt, units, react_ext));
    }
    return *this;
}

auto Material::elements() const -> ElementMap
{
    ElementMap total;
    for (const auto& c : constituents_) {
        for (const auto& kv : c.elementMoles()) {
            total[kv.first] += kv.second;
        }
    }
    return total;
}

auto Material::b() const -> Vector
{
    if (!adapter_) {
        throw std::runtime_error(
            "xGEMS::Material::b(): no engine bound. "
            "Construct with Material(engine).");
    }
    const auto& names = adapter_->elementNames();
    Vector bvec = Vector::Constant(static_cast<Index>(names.size()), default_amount_);

    for (const auto& kv : elements()) {
        const Index i = adapter_->indexOf(kv.first);
        if (i < 0) {
            throw std::runtime_error(
                "xGEMS::Material::b(): element '" + kv.first +
                "' from recipe is not in the engine's element list.");
        }
        bvec[i] += kv.second;
    }
    return bvec;
}

auto Material::bMap() const -> ElementMap
{
    if (!adapter_) {
        throw std::runtime_error(
            "xGEMS::Material::bMap(): no engine bound.");
    }
    ElementMap result;
    for (const auto& nm : adapter_->elementNames()) result[nm] = default_amount_;

    for (const auto& kv : elements()) {
        auto it = result.find(kv.first);
        if (it == result.end()) {
            throw std::runtime_error(
                "xGEMS::Material::bMap(): element '" + kv.first +
                "' from recipe is not in the engine's element list.");
        }
        it->second += kv.second;
    }
    return result;
}

auto Material::massKg() const -> double
{
    double m = 0.0;
    for (const auto& c : constituents_) m += c.massKg();
    return m;
}

auto Material::moles() const -> double
{
    double n = 0.0;
    for (const auto& c : constituents_) n += c.moles();
    return n;
}

auto Material::total(const std::string& units) const -> double
{
    checkUnit(units);
    if (auto it = MOLE_TO_MOL.find(units); it != MOLE_TO_MOL.end()) {
        return moles() / it->second;
    }
    return massKg() / MASS_TO_KG.at(units);
}

auto Material::scale(double factor) -> Material&
{
    for (auto& c : constituents_) c.value *= factor;
    return *this;
}

auto Material::scaleToMass(double kg) -> Material&
{
    const double current = massKg();
    if (current <= 0.0) {
        throw std::runtime_error(
            "xGEMS::Material::scaleToMass: current mass is zero.");
    }
    return scale(kg / current);
}

auto Material::scaleToMoles(double mol) -> Material&
{
    const double current = moles();
    if (current <= 0.0) {
        throw std::runtime_error(
            "xGEMS::Material::scaleToMoles: current total moles is zero.");
    }
    return scale(mol / current);
}

auto Material::setQuantity(double value, const std::string& units) -> Material&
{
    const double current = total(units);
    if (current <= 0.0) {
        throw std::runtime_error(
            "xGEMS::Material::setQuantity: cannot rescale a recipe whose "
            "current total is non-positive.");
    }
    return scale(value / current);
}

auto Material::operator()(double value, const std::string& units) const -> Material
{
    const double current = total(units);
    if (current <= 0.0) {
        throw std::runtime_error(
            "xGEMS::Material::operator(): cannot scale a recipe whose "
            "current total is non-positive.");
    }
    Material copy(*this);
    std::ostringstream os;
    os << name_ << "(" << value << " " << units << ")";
    copy.name_ = os.str();
    copy.scale(value / current);
    return copy;
}

auto Material::operator+(const Material& other) const -> Material
{
    Material result(*this);
    if (!result.adapter_)
        result.adapter_ = other.adapter_;
    result.name_    = "(" + name_ + " + " + other.name_ + ")";
    result.constituents_.insert(
        result.constituents_.end(),
        other.constituents_.begin(), other.constituents_.end());
    return result;
}

auto Material::toString() const -> std::string
{
    std::ostringstream os;
    os << "Material('" << name_ << "', "
       << constituents_.size() << " constituents, "
       << (massKg() * 1.0e3) << " g, "
       << moles() << " mol)";
    return os.str();
}

} // namespace xGEMS
