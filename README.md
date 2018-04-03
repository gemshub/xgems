
## xGEMS 

This is a numerical solver of chemical equilibria for thermodynamic modelling, extended (relative to GEMS3K code) with new C++ and Python APIs (compatible with [Reaktoro framework](http://reaktoro.org)), and supposed to replace GEMS3K in next-generation software.  

### Briefly about xGEMS

* Implements the improved GEM IPM-3 algorithm with excellent mass balance precision and fast convergence to Gibbs energy minimum even in very complex non-ideal chemical systems with two-sided metastability constraints.
* Written in C/C++. Using compiler directives, the xGEMS code can be compiled as a standalone program e.g. 'gemcalc'; as a static or dynamic library for coupling with the mass transport simulator or another code; or as part of the GEM-Selektor v4 code together with GUI and databases.
* Version: currently 4.0.0.
* Will be distributed under the terms of LGPL v.3 license. 

### How to download xGEMS source code

* In your home directory, make a folder named e.g. ~/git/xGEMS.
* Change into ~/git/xGEMS and clone this repository from https://bitbucket.org/gems4/xgems.git using a preinstalled free git client, e.g. SourceTree. 
* Alternatively on Mac OS X or linux, open a terminal and type in the command line:
~~~
cd ~/git/xGEMS
git clone https://bitbucket.org/gems4/xgems.git .
~~~


### How to build xGEMS library and examples

* You will need to install Eigen [(the development branch)](http://bitbucket.org/eigen/eigen/get/default.tar.bz2). 
* We assume you have unpacked it to your home directory and renamed to ~/unpacked-eigen-dir.
~~~
cd ~/unpacked-eigen-dir
mkdir build 
cd build
cmake ..
sudo make install
~~~

* You will also need to install Pybind11 (v2.2.2):
~~~
cd Downloads
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build
cd build
cmake .. -DPYBIND11_TEST=OFF
sudo make install
~~~

* To compile xGEMS and demos:
~~~
cd ~/git/xGEMS
mkdir build
cd build
cmake .. -DPYTHON_EXECUTABLE=/usr/bin/python2.7
make -j 3
make demos
~~~
If you use Python 3.5, in the anove commands, change "python2.7" to "python3.5". 
For compiling in the debug mode, add -DCMAKE_BUILD_TYPE=Debug as cmake parameter.

* To execute the c++ demo:
~~~
cd ~/git/xGEMS/build/bin
./demo1
~~~

* To execute the python demo:
* First, make sure the xGEMS Python module can be found by Python. 
* This can be done by setting the `PYTHONPATH` environment variable to the path where the xGEMS module is located:

~~~
export PYTHONPATH=/home/username/pathto-xGEMS/build/lib
~~~

* Then run the demo - if built with python 2.7:
~~~
 cd ~git/xGEMS/demos/
 python demo1.py
~~~

* If built with python 3.5: in demo1.py, change the last line from "print chemicalengine" to "print(chemicalengine)", save and execute:
~~~
 cd ~git/xGEMS/demos/
 python3 demo1.py
~~~



There are yet things to do.


## Edit a file

You’ll start by editing this README file to learn how to edit a file in Bitbucket.

1. Click **Source** on the left side.
2. Click the README.md link from the list of files.
3. Click the **Edit** button.
4. Delete the following text: *Delete this line to make a change to the README from Bitbucket.*
5. After making your change, click **Commit** and then **Commit** again in the dialog. The commit page will open and you’ll see the change you just made.
6. Go back to the **Source** page.

---

**Edit a file, create a new file, and clone from Bitbucket in under 2 minutes**

When you're done, you can delete the content in this README and update the file with details for others getting started with your repository.

*We recommend that you open this README in another tab as you perform the tasks below. You can [watch our video](https://youtu.be/0ocf7u76WSo) for a full demo of all the steps in this tutorial. Open the video in a new tab to avoid leaving Bitbucket.*

---


## Create a file

Next, you’ll add a new file to this repository.

1. Click the **New file** button at the top of the **Source** page.
2. Give the file a filename of **contributors.txt**.
3. Enter your name in the empty file space.
4. Click **Commit** and then **Commit** again in the dialog.
5. Go back to the **Source** page.

Before you move on, go ahead and explore the repository. You've already seen the **Source** page, but check out the **Commits**, **Branches**, and **Settings** pages.

---

## Clone a repository

Use these steps to clone from SourceTree, our client for using the repository command-line free. Cloning allows you to work on your files locally. If you don't yet have SourceTree, [download and install first](https://www.sourcetreeapp.com/). If you prefer to clone from the command line, see [Clone a repository](https://confluence.atlassian.com/x/4whODQ).

1. You’ll see the clone button under the **Source** heading. Click that button.
2. Now click **Check out in SourceTree**. You may need to create a SourceTree account or log in.
3. When you see the **Clone New** dialog in SourceTree, update the destination path and name if you’d like to and then click **Clone**.
4. Open the directory you just created to see your repository’s files.

Now that you're more familiar with your Bitbucket repository, go ahead and add a new file locally. You can [push your change back to Bitbucket with SourceTree](https://confluence.atlassian.com/x/iqyBMg), or you can [add, commit,](https://confluence.atlassian.com/x/8QhODQ) and [push from the command line](https://confluence.atlassian.com/x/NQ0zDQ).