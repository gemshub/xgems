
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

* Install Dependencies by executing in ```~/xGEMS$```

```sh
#!bash
sudo ./install-dependencies.sh
```

* To build xGEMS and install it in your home directory or in the system directory (as in the example below), a typical sequence of commands can be executed in the terminal:

```sh
cd ~/xGEMS
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
```

* To execute the c++ demo:

```sh
cd ~/build/bin
./demo1
```

* To execute the python demo: 

```sh
 cd ~git/xGEMS/demos/
 python3 demo1.py
 ```

There are yet things to do.
