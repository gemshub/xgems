# xGEMS

[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/xgems?style=for-the-badge&logo=conda-forge)](https://anaconda.org/conda-forge/xgems)
[![Docs](https://img.shields.io/readthedocs/xgems?style=for-the-badge)](https://xgems.readthedocs.io)
[![Binder](https://img.shields.io/badge/try%20it-binder-579ACA?style=for-the-badge)](https://mybinder.org/v2/gh/gemshub/xgems-jupyter/main?urlpath=lab/tree/how-to-use-xgems-examples.ipynb)

xGEMS calculates chemical equilibrium in complex non-ideal systems by Gibbs energy minimization.
It gives the [GEMS3K](https://github.com/gemshub/GEMS3K) solver a C++ and Python interface.
It runs on Linux, macOS and Windows.

## Example

Calcite in 1 kg of water. More CO2 makes the water more acidic, and the calcite dissolves:

```python
import xgems

engine = xgems.ChemicalEngineDicts("demos/resources/CalciteBC/CalciteBC-dat.lst")
engine.T, engine.P = 298.15, 1e5            # K, Pa

for co2 in [0.0, 0.001, 0.01, 0.1]:         # mol of CO2 added
    engine.set_bulk_composition({"H": 111.017, "O": 55.508 + 0.03 + 2*co2,
                                 "Ca": 0.01, "C": 0.01 + co2})
    engine.equilibrate()
    print(f"CO2 {co2:5.3f} mol | pH {engine.pH:5.2f} | "
          f"calcite left {engine.phases_moles['Calcite']:.4f} mol | "
          f"Ca in water {engine.aq_elements_molality['Ca']:.2e} mol/kg")

# Optima solver: set the pH to 5.5 and get the acid that it needs
engine.setSolverMode("aop")
engine.setpHTarget(5.5)
engine.equilibrate()
print(f"pH {engine.pH:.2f} needs {-engine.controlConditionTitrant('pH'):.2e} mol of H+")
```

```text
CO2 0.000 mol | pH 10.31 | calcite left 0.0095 mol | Ca in water 3.12e-04 mol/kg
CO2 0.001 mol | pH  7.89 | calcite left 0.0092 mol | Ca in water 8.49e-04 mol/kg
CO2 0.010 mol | pH  6.53 | calcite left 0.0053 mol | Ca in water 4.67e-03 mol/kg
CO2 0.100 mol | pH  6.01 | calcite left 0.0004 mol | Ca in water 9.64e-03 mol/kg
pH 5.50 needs 1.40e-02 mol of H+
```

The chemical system (species, phases, thermodynamic data) comes from a project that
[GEMS-Selektor](https://gems.web.psi.ch) exports.

More examples as Jupyter notebooks: [xgems-jupyter](https://github.com/gemshub/xgems-jupyter).
[Try them in your browser with Binder](https://mybinder.org/v2/gh/gemshub/xgems-jupyter/main?urlpath=lab/tree/how-to-use-xgems-examples.ipynb),
no installation needed.

## Features

* Two solvers. The GEMS3K IPM solver is fast and is the default. The Optima solver is more reliable for solid
  solutions, melts and miscibility gaps. Select the solver for each calculation with `setSolverMode()`.
* Cold start, or warm start from the previous result. A warm start costs much less in sequences of small changes,
  for example in reactive transport.
* Target pH and Eh: the solver finds the acid, base or electrons that give the target.
* A failed calculation is tried again in steps before xGEMS reports it.
* Results by name: amounts, activities, saturation indices, volumes, masses and thermodynamic properties of
  species, phases and the system.
* Use from Python or C++, or link the library into a transport code.

## Solver modes

* `setSolverMode()` selects the solver for each calculation: `"aia"`, `"sia"` (GEMS3K IPM solver, cold and warm
  start), `"aop"`, `"sop"` (Optima solver, cold and warm start), `"hop"`, `"shp"` (GEMS3K IPM solver, then an Optima
  check).
* `solverMode()` tells which mode runs. `builtWithOptima()` tells if the Optima modes are available.
* If `warmstart` is set, `"aia"` runs as `"sia"`, `"aop"` as `"sop"` and `"hop"` as `"shp"`.
* The default is `"aia"`.

## pH and Eh targets

* `setpHTarget()` and `setEhTarget()` set a target pH or Eh in the Optima modes. The solver adds the needed amount
  of acid, base or electrons.
* `controlConditionTitrant()` gives that amount. For a pH target it is the amount of H+ removed, so a negative value
  means H+ added. `clearControlConditions()` removes the targets.

## Documentation

The full documentation is on [xgems.readthedocs.io](https://xgems.readthedocs.io): installation, solver guide,
settings, advanced usage, and the Python and C++ API. More examples are in [demos/](demos) and in the
[Jupyter notebooks](https://github.com/gemshub/xgems-jupyter).

## Installation with conda

```bash
conda install -c conda-forge xgems
```

To get conda, install [Miniforge](https://conda-forge.org/download/).

## Build from source

Get the source code:

```bash
git clone https://github.com/gemshub/xgems.git
cd xgems
```

### Option 1: in a conda environment (recommended)

The environment file installs all dependencies, including GEMS3K with the Optima solver.

```bash
conda install -c conda-forge conda-devenv
conda devenv
conda activate xgems
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
make
make install
```

### Option 2: without conda

The script builds the dependencies and GEMS3K. It does not build Optima, so the Optima modes use the GEMS3K IPM
solver.

```bash
sudo ./install-dependencies.sh
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
```

### Check the installation

```bash
cd build/bin
./demo1             # C++ demo
python3 demo1.py    # Python demo
```

## License

xGEMS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version. See [COPYING.LESSER](COPYING.LESSER) and [COPYING](COPYING).

xGEMS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

xGEMS uses GEMS3K (LGPL-3.0-or-later) and the modified Optima of gemshub (LGPL-2.1-or-later, see its LICENSE.md).

## For developers: build the documentation

First build and install xGEMS. Then:

```bash
conda install -c conda-forge doxygen sphinx breathe sphinx_rtd_theme
doxygen Doxyfile
sphinx-build -b html docs/sphinx/source docs/build/html
```
