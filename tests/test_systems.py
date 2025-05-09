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
import unittest
from xgems import *
from numpy import *

def gems3k_equilibrate(fullpath, T=None, P=None):
    engine = ChemicalEngine(fullpath)

    T1 = engine.temperature() if T is None else T
    P1 = engine.pressure() if P is None else P
    b = engine.elementAmounts()

    ret = engine.equilibrate(T1, P1, b)
    #print("Run: ", fullpath, " result: ", ret, "time: ", engine.elapsedTime())
    return ret

def gems3k_folder(fullpath):
    for x in os.listdir(fullpath):
        if x.endswith("-dat.lst"):
            return os.path.join(fullpath, x)
    return None

class TestCalculations(unittest.TestCase):

    def setUp(self):
        path = "gems3k"
        self.dirs = [os.path.join(path, name) for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
        print(self.dirs)


    def test_equilibrate(self):
        for folder in self.dirs:
            with self.subTest(msg = folder):
                lst_file = gems3k_folder(folder)
                if lst_file is not None:
                    self.assertEqual(gems3k_equilibrate(lst_file), 3)

    def test_equilibrate_1_25(self):
        for folder in self.dirs:
            with self.subTest(msg = folder):
                lst_file = gems3k_folder(folder)
                if lst_file is not None:
                    self.assertEqual(gems3k_equilibrate(lst_file, 298.15, 100000), 2)


if __name__ == '__main__':
    unittest.main()


