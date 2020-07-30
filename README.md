
## xGEMS 

This is a numerical solver of chemical equilibria for thermodynamic modelling, extended (relative to GEMS3K code) with new C++ and Python APIs (compatible with [Reaktoro framework](http://reaktoro.org)), and supposed to replace GEMS3K in next-generation software.  

### Briefly about xGEMS

* Currently uses the [GEMS3K](https://bitbucket.org/gems4/gems3k) code that implements the improved GEM IPM-3 algorithm for Gibbs energy minization in very complex non-ideal chemical systems with two-sided metastability constraints.
* Written in C/C++. The xGEMS code can be compiled as a standalone program (see /demos); or built as a static or dynamic library for coupling with a mass transport simulator or another code.

### License
* GEMS3K is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* GEMS3K is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
* You should have received a copy of the GNU Lesser General Public License along with GEMS3K code. If not, see http://www.gnu.org/licenses/. 

### How to download xGEMS source code

* In your home directory, make a folder named e.g. ~/git/xGEMS.
* Change into ~/git/xGEMS and clone this repository from https://bitbucket.org/gems4/xgems.git using a preinstalled free git client, e.g. SourceTree. 
* Alternatively on Mac OS X or linux, open a terminal and type in the command line:
~~~
cd ~/git/xGEMS
git clone https://bitbucket.org/gems4/xgems.git .
~~~

### How to build xGEMS library and examples

xGEMS uses GEMS3K and Egien3 as external projects. These are downloaded and compiled together with xGEMS during the make step. 

For compiling xGEMS python bindings you need a recent version of pybind11 installed (including development headers, 2.2.0 is the minimum version required). In case your distribution does not provide pybind11 development packages, or the version is to old, you can install it as follows:

~~~
mkdir -p ~/code
cd ~/code
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build
cd build
cmake .. -DPYBIND11_TEST=OFF -DCMAKE_INSTALL_PREFIX=/your/local/codes/directory
make install
~~~

CMake uses `find_package` to search for pybind11 `pybind11Config.cmake` configuration file. Installing pybind11 using `pip install` or other ways might not provide the required configuration file and will result in errors.

* To compile xGEMS with python bindings and demos:

`-DCMAKE_PREFIX_PATH=/your/local/codes/directory` is necessary if pybind11 is not installed in the default include path `/usr/include/` or `/usr/lib/cmake/`, provided that `/your/local/codes/directory/pybind11/pybind11Config.cmake` cmake config file is found at this path.

~~~
cd ~/git/xGEMS
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=/your/local/codes/directory
make -j 3
make demos
~~~

* For building xGEMS for use from C++ codes only (i.e. no python in the system), replace the above cmake command with
~~~
cmake .. -DXGEMS_BUILD_PYTHON=OFF
~~~

* For pointing to a specific python version use `-DPYTHON_EXECUTABLE=/usr/bin/python3.6` or `-DPYTHON_EXECUTABLE=$(which python)` to get the default python path.

* For compiling in the debug mode, add -DCMAKE_BUILD_TYPE=Debug as cmake parameter.

* To install xGEMS into /usr/local as a library with /includes:
~~~
cd ~/git/xGEMS/build
sudo make install 
~~~

* To execute the c++ demo:
~~~
cd ~/git/xGEMS/build/bin
./demo1
~~~

* To execute the python demo: first, make sure the xGEMS Python module can be found by Python. 
* This can be done by setting the `PYTHONPATH` environment variable to the path where the xGEMS module is located:
~~~
export PYTHONPATH=/home/username/pathto-xGEMS/build/lib
~~~

(use your username instead of "username" and actual path to xGEMS e.g. git/xGEMS instead of "pathto-xGEMS"   

* Then run the demo - if built with python 2.7:
~~~
 cd ~git/xGEMS/demos/
 python demo1.py
~~~

* If built with python 3.6: in demo1.py, change the last line from "print chemicalengine" to "print(chemicalengine)", save and execute:
~~~
 cd ~git/xGEMS/demos/
 python3 demo1.py
~~~

There are yet things to do.
