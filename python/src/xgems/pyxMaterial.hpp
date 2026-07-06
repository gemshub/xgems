// xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2018-2026 G.D. Miron, D. Kulik, S. Dmytriieva and contributors
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

#pragma once

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

// xGEMS includes
#include <xGEMS/Material.hpp>
#include <xGEMS/ChemicalEngine.hpp>
#include <xGEMS/ChemicalEngineMaps.hpp>
using namespace xGEMS;


inline void exportMaterial(py::module& m)
{
    // ─── Free helper ────────────────────────────────────────────────
    m.def("parse_chemical_formula", &parseChemicalFormula,
          py::arg("formula"),
          R"doc(
Parse a chemical-formula string into a {element: count} dictionary.

Supports nested parentheses and hydrate notation with '.' separators.

:param str formula: Chemical formula (e.g. "CaSO4.2H2O", "Ca(OH)2", "Cl").
:return dict: {element_symbol: float count per formula unit}

**Example:**

.. code-block:: python

    parse_chemical_formula("CaSO4.2H2O")
    # → {"Ca": 1, "S": 1, "O": 6, "H": 4}
)doc");

    // ─── Material class ────────────────────────────────────────────
    py::class_<Material>(m, "Material",
        R"doc(
A recipe of constituents bound to a ChemicalEngine or ChemicalEngineDicts.

Build with .add(), scale by calling ``brine(1.0, "kg")``, rescale in
place with ``brine.set_quantity(1.0, "kg")``, and combine recipes with
``+``. The result of ``.b()`` is an Eigen vector indexed by the
engine's element list — ready to pass to ``equilibrate``.

Supported units:
  - Mass : "kg", "g", "mg", "ug"
  - Moles: "mol" (or "moles"), "mmol", "umol"

**Example:**

.. code-block:: python

    from xgems import ChemicalEngine, Material

    engine = ChemicalEngine()
    engine.initialize("system-dat.lst")

    # Seawater (Millero et al.) — H2O separated from dissolved ions
    seawater = Material(engine, "seawater")
    seawater.add("H2O", 1.00, "kg")
    seawater.add(
        {"Cl": 0.5467499,  "Na": 0.4689999,
         "Mg": 0.05282088, "Ca": 0.01037339,
         "K":  0.01020831, "S":  0.02823655,
         "C":  0.001966361},
        1.0, "mol",
    )

    # Basalt glass — oxide blend
    basalt = Material(engine, "basalt_glass")
    for ox, mol in [("SiO2", 1.000), ("Al2O3", 0.175), ("FeO",   0.169),
                    ("CaO",  0.270), ("MgO",   0.260), ("Na2O",  0.040),
                    ("K2O",  0.004), ("TiO2",  0.020), ("Fe2O3", 0.006)]:
        basalt.add(ox, mol, "mol")

    # Leucogranite — mineral assemblage by mass
    granite = Material(engine, "leucogranite")
    granite.add("KAlSi3O8",        17.0, "g")    # Microcline
    granite.add("KAl3Si3O10(OH)2", 19.0, "g")    # Muscovite
    granite.add("NaAlSi3O8",       29.0, "g")    # Albite
    granite.add("SiO2",            35.0, "g")    # Quartz

    # Combine and scale
    mix = seawater(1.0, "kg") + basalt(10.0, "g")
    engine.equilibrate(298.15, 1.0e5, mix.b())

    # In-place rescale (vs the copy returned by ``granite(0.5, "kg")``)
    granite.set_quantity(0.5, "kg")
)doc")

        // Constructors -------------------------------------------------
        .def(py::init<>(),
             R"doc(Default constructor; not yet bound to an engine.

``.b()`` will throw until the Material is given an engine (e.g. by
being the result of ``+`` with a bound Material).
)doc")

        .def(py::init<const ChemicalEngine&, std::string>(),
             py::arg("engine"), py::arg("name") = "material",
             py::keep_alive<1, 2>(),
             R"doc(
Construct an empty Material bound to a ChemicalEngine.

:param ChemicalEngine engine: Engine whose element list and molar masses
    are used. The engine must outlive the Material (pybind keeps the
    Python-side reference alive automatically via ``keep_alive``).
:param str name: Optional display name.
)doc")

        .def(py::init<const ChemicalEngineMaps&, std::string>(),
             py::arg("engine"), py::arg("name") = "material",
             py::keep_alive<1, 2>(),
             R"doc(
Construct an empty Material bound to a ChemicalEngineDicts.

Same as the ChemicalEngine overload but uses the dict-style API. The
Python ``ChemicalEngineDicts`` class is the pybind-side name for the
C++ ``ChemicalEngineMaps`` type.

:param ChemicalEngineDicts engine: Engine to bind to.
:param str name: Optional display name.
)doc")

        // .add() overloads ---------------------------------------------
        .def("add",
             py::overload_cast<const std::string&, double,
                               const std::string&, double>(&Material::add),
             py::arg("formula"), py::arg("value") = 1.0,
             py::arg("units") = "mol", py::arg("react_ext") = 1.0,
             py::return_value_policy::reference_internal,
             R"doc(
Add a constituent by chemical-formula string.

:param str formula: e.g. "CaCO3", "H2O", "Cl", "CaSO4.2H2O".
:param float value: Amount in ``units`` (default 1.0).
:param str units: One of "mol"/"moles", "mmol", "umol", "kg", "g", "mg", "ug".
:param float react_ext: Reaction-extent multiplier in [0,1] (default 1).
:return Material: ``self`` for chaining.
)doc")

        .def("add",
             [](Material& self,
                const ElementMap& formula,
                py::object value_obj,
                const std::string& units,
                double react_ext) -> Material& {
                 if (value_obj.is_none()) {
                     return self.add(formula, units, react_ext);
                 }
                 // unit string passed as second positional arg (absolute-amount mode)
                 if (py::isinstance<py::str>(value_obj)) {
                     return self.add(formula, value_obj.cast<std::string>(), react_ext);
                 }
                 return self.add(formula, value_obj.cast<double>(), units, react_ext);
             },
             py::arg("formula"),
             py::arg("value") = py::none(),
             py::arg("units") = "mol",
             py::arg("react_ext") = 1.0,
             py::return_value_policy::reference_internal,
             R"doc(
Add a constituent by an explicit element-count dict.

When ``value`` is omitted (or ``None``), the dict values are treated as
the **absolute amounts** of each element in ``units`` — each entry
becomes its own single-element constituent and the total is the sum of
all values after unit conversion.

When ``value`` is a float, the dict values are element counts per
formula unit and the whole group is scaled by ``value`` in ``units``
(the original behaviour).

:param dict formula: element → amount or element → count per formula unit.
:param value: Scale factor in ``units``. Omit (or pass ``None``) to use
    absolute-amount mode.
:type value: float or None
:param str units: One of "mol"/"moles", "mmol", "umol", "kg", "g", "mg", "ug".
:param float react_ext: Reaction-extent multiplier in [0,1] (default 1).
:return Material: ``self`` for chaining.

**Absolute-amount mode** (value omitted):

.. code-block:: python

    sw.add("H2O", 1.00, "kg")
    sw.add(
        {"Cl": 0.5467499, "Na": 0.4689999,
         "Mg": 0.05282088, "Ca": 0.01037339,
         "K":  0.01020831, "S":  0.02823655,
         "C":  0.001966361},
        "mol",               # each value IS the amount in mol
    )

**Scale mode** (value provided):

.. code-block:: python

    sw.add({"Na": 1, "Cl": 1}, 1.0, "mol")   # 1 mol of NaCl
)doc")

        // Scaling syntax: material(value, units) returns a copy --------
        .def("__call__", &Material::operator(),
             py::arg("value") = 1.0, py::arg("units") = "mol",
             R"doc(
Return a scaled COPY whose total equals ``(value, units)``.

The original Material is unchanged. Useful for the natural-looking
``brine(1.0, "kg") + rock(10.0, "g")`` syntax.
)doc")

        .def("__add__", &Material::operator+,
             R"doc(Concatenate two recipes into a new Material.)doc")

        // In-place scaling ---------------------------------------------
        .def("scale", &Material::scale, py::arg("factor"),
             py::return_value_policy::reference_internal,
             R"doc(Multiply every constituent's value by ``factor`` in place.)doc")

        .def("scale_to_mass", &Material::scaleToMass, py::arg("kg"),
             py::return_value_policy::reference_internal,
             R"doc(Rescale in place so the total mass equals ``kg`` kilograms.)doc")

        .def("scale_to_moles", &Material::scaleToMoles, py::arg("mol"),
             py::return_value_policy::reference_internal,
             R"doc(Rescale in place so the total amount equals ``mol`` moles.)doc")

        .def("set_quantity", &Material::setQuantity,
             py::arg("value"), py::arg("units") = "mol",
             py::return_value_policy::reference_internal,
             R"doc(
Rescale in place so the total equals ``(value, units)``.

Companion to ``material(value, units)`` which returns a scaled copy.
``set_quantity`` modifies the Material directly and returns ``self``
for chaining.
)doc")

        // Aggregates ---------------------------------------------------
        .def("elements", &Material::elements,
             R"doc(
Return the recipe's element-mole contribution as a dict.

:return dict[str,float]: {element_name: total moles in this material}.
)doc")

        .def("b", &Material::b,
             R"doc(
Return the bulk-composition vector ``b`` indexed by the engine's
element list (length ``engine.numElements()``).

Each non-``Zz`` engine element is clamped to at least ``min_amount``
(default 1e-15). The charge element ``Zz`` is always 0. Ready to be passed to
``engine.equilibrate(T, P, b)``.

**Example:**

.. code-block:: python

    engine.equilibrate(298.15, 1e5, material.b())
)doc")

        .def("b_dict", &Material::bMap,
             R"doc(
Return the bulk-composition as a ``{element_name: moles}`` dict.

Natural form for ChemicalEngineDicts users. Includes one entry per
engine element; each non-``Zz`` value is clamped to at least
``min_amount`` (default 1e-15). ``Zz`` is always 0.

**Example:**

.. code-block:: python

    engine.equilibrate(298.15, 1e5, material.b_dict())
)doc")

        .def("mass_kg", &Material::massKg,
             R"doc(Total mass of the recipe in kg.)doc")
        .def("moles",   &Material::moles,
             R"doc(Total amount of the recipe in mol.)doc")
        .def("total",   &Material::total, py::arg("units") = "mol",
             R"doc(Total in any supported unit string.)doc")

        // Inspection ---------------------------------------------------
        .def_property("name", &Material::name,
                      [](Material& self, const std::string& n) { self.setName(n); },
                      "Display name of the Material.")
        .def_property("min_amount",
                      &Material::minAmount,
                      [](Material& self, double v) { self.setMinAmount(v); },
                      R"doc(
Minimum [mol] returned for each non-``Zz`` entry in ``b()`` and ``b_dict()``.

The GEM solver requires strictly positive bulk amounts. Any non-``Zz``
entry below this value is raised to the floor (default 1e-15). Set to
0.0 to disable the clamp and pass exact recipe amounts to the solver.

**Example:**

.. code-block:: python

    material.min_amount = 1e-12   # coarser floor
    material.min_amount = 0.0     # disable floor
    print(material.min_amount)    # 0.0
)doc")
        .def("__len__", &Material::size)
        .def("__repr__", &Material::toString)
        ;
}
