
## xGEMS 

Linux, OSX, Windows

This is a numerical solver of chemical equilibria for thermodynamic modelling, extended (relative to GEMS3K code) with new [C++ and Python APIs](https://gemshub.github.io/site/start/gemstandalone/documentation/#cpython-chemicalengine-interface), and supposed to replace GEMS3K in next-generation software.  

### Documentation

For detailed information about how to use this project through its C++ and Python interface interface [View Documentation](https://xgems.readthedocs.io)

### Briefly about xGEMS

* Currently uses the [GEMS3K](https://github.com/gemshub/GEMS3K) code that implements the improved GEM IPM-3 algorithm for Gibbs energy minimization in very complex non-ideal chemical systems with two-sided metastability constraints.
* Written in C/C++. The xGEMS code can be compiled as a standalone program (see /demos); or built as a static or dynamic library for coupling with a mass transport simulator or another code.

### License
* GEMS3K is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* GEMS3K is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
* You should have received a copy of the GNU Lesser General Public License along with GEMS3K code. If not, see http://www.gnu.org/licenses/. 

### Installation using Conda 

[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/xgems?style=for-the-badge&logo=conda-forge)](https://anaconda.org/conda-forge/xgems)

xGEMS can be easily installed using [Conda](https://www.anaconda.com/docs/getting-started/miniconda/install) package manager. If you have Conda or Miniforge installed, first add the conda-forge channel by executing 

```bash
conda config --add channels conda-forge
```

install xGEMS by executing the following command:

```bash
conda install xgems
```

Conda can be installed from [Miniconda](https://conda.io/miniconda.html).

### How to download xGEMS source code

* In your home directory, make a folder named e.g. ~/git/xGEMS.
* Change into ~/git/xGEMS and clone this repository from https://bitbucket.org/gems4/xgems.git using a preinstalled free git client, e.g. SourceTree. 
* Alternatively on Mac OS X or linux, open a terminal and type in the command line:
```bash
cd ~/git/xGEMS
git clone https://bitbucket.org/gems4/xgems.git .
```

### How to build xGEMS library and examples

#### **Option 1: Build xGEMS in a Conda Environment Using `conda-devenv`**

1. **Install Conda and `conda-devenv`**  
   Ensure you have Conda installed. If not, download and install Miniconda or Anaconda from [here](https://docs.conda.io/en/latest/miniconda.html).  
   Install `conda-devenv`:
   ```bash
   conda install -c conda-forge conda-devenv
   ```

2. **Create and Activate the Conda Environment** 
   ```bash
   conda devenv
   conda activate xgems
   ```
3. **Build xGEMS**
   ```bash
   cd ~/git/xGEMS
   mkdir build
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
   make
   make install
   ```

#### **Option 2: Build xGEMS without Conda**


1. **Install Dependencies**

   by executing in ```~/xGEMS$```
   
   ```bash
   sudo ./install-dependencies.sh
   ```

2. **Build xGEMS** 

   To build xGEMS and install it in your home directory or in the system directory (as in the example below), a typical sequence of commands can be executed in the terminal:
   
   ```sh
   cd ~/xGEMS
   mkdir build
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
   make
   sudo make install
   ```

#### Verify installation

   * To execute the c++ demo:
   
   ```bash
   cd ~/build/bin
   ./demo1
   ```
   
   * To execute the python demo: 
   
   ```bash
    cd ~git/xGEMS/demos/
    python3 demo1.py
    ```

## For developers

### Update documentation 

Install doxygen

```bash
sudo apt install doxygen  # Linux
```
Install addons 
```bash
pip install sphinx breathe sphinx_rtd_theme
pip install sphinx-autodoc-typehints
pip install sphinxcontrib-autodoc-doxygen
```

(first build and install xgems, then update documentation)

Generate xml and html
```bash
doxygen Doxyfile
sphinx-build -b html docs/sphinx/source docs/html
```




