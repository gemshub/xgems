# xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
#
# Copyright (C) 2018 Allan Leal, Dmitrii Kulik
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import os
from xgems import *
from numpy import *



def gems3k_equilibrate(fullpath, T=None, P=None):
    try:
        engine = ChemicalEngine(fullpath)

        T1 = engine.temperature() if T is None else T
        P1 = engine.pressure() if P is None else P
        b = engine.elementAmounts()

        ret = engine.equilibrate(T1, P1, b)
        print("Run: ", fullpath, " result: ", ret, "time: ", engine.elapsedTime())
    except RuntimeError:
        print("Illegal input: ", fullpath)

def gems3k_folder(name, fullpath, T=None, P=None):
    #print(name + " " + fullpath)
    for x in os.listdir(fullpath):
        if x.endswith("-dat.lst"):
            gems3k_equilibrate(os.path.join(fullpath, x), T, P)


def gems3k_folders(path=".", T=None, P=None):
    for name in os.listdir(path):
        thepath = os.path.join(path, name)
        if os.path.isdir(thepath):
            gems3k_folder(name, thepath, T, P)


if __name__ == "__main__":
    #connect all loggers
    #update_loggers(True, "test_demo1.log", 0)
    gems3k_folders("gems3k")



