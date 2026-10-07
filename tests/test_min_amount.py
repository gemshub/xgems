# xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
#
# Checks the default minimum amount given to elements that the bulk composition does not contain (1e-11 mol).
# Below the Optima solver's floor, GEMS3K warns "N element(s) have less material than the solver's floor amount
# allows" and repairs the mass balance afterwards; the default is chosen above that for ordinary aqueous systems.

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

from xgems import ChemicalEngine, ChemicalEngineDicts, Material

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.join(HERE, "gems3k", "CASHNK1", "1TCi+CH-dat.lst")

OPTIMA = ChemicalEngine.builtWithOptima()


def recipe(engine):
    """1 kg of water with a little Ca and Si: the element Nit (nitrogen) is not in the recipe."""
    material = Material(engine, "m")
    material.add("H2O", 1.0, "kg")
    material.add({"Ca": 1e-3, "Si": 1e-3, "O": 3e-3})
    return material


class TestDefaultMinAmount(unittest.TestCase):

    def test_material_default(self):
        engine = ChemicalEngineDicts(PROJECT)
        material = recipe(engine)
        self.assertEqual(material.min_amount, 1e-11)
        self.assertEqual(material.b_dict()["Nit"], 1e-11)      # the absent element is at the default

    def test_material_value_can_be_changed(self):
        engine = ChemicalEngineDicts(PROJECT)
        material = recipe(engine)
        material.min_amount = 1e-15
        self.assertEqual(material.b_dict()["Nit"], 1e-15)
        material.min_amount = 0.0                                # no floor at all
        self.assertNotIn("Nit", {k: v for k, v in material.b_dict().items() if v > 0})

    def test_engine_default(self):
        engine = ChemicalEngineDicts(PROJECT)
        engine.clear()                                           # every element at the minimum amount
        self.assertEqual(engine.bulk_composition["Nit"], 1e-11)
        self.assertEqual(engine.bulk_composition["Ca"], 1e-11)


@unittest.skipUnless(OPTIMA, "GEMS3K was built without Optima")
class TestFloorWarning(unittest.TestCase):

    CODE = """
    import os, sys, xgems
    from xgems import ChemicalEngineDicts, Material
    engine = ChemicalEngineDicts({project!r})
    xgems.update_loggers(False, "log.txt", 3)
    engine.setSolverMode("aop")
    material = Material(engine, "m")
    material.add("H2O", 1.0, "kg")
    material.add({{"Ca": 1e-3, "Si": 1e-3, "O": 3e-3}})
    {setting}
    engine.equilibrate(298.15, 1e5, material)
    """

    def floor_warnings(self, setting):
        with tempfile.TemporaryDirectory() as cwd:
            code = textwrap.dedent(self.CODE).format(project=PROJECT, setting=setting)
            proc = subprocess.run([sys.executable, "-c", code], cwd=cwd, capture_output=True, text=True, timeout=300)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            path = os.path.join(cwd, "log.txt")
            return open(path).read().count("element(s) have less material") if os.path.exists(path) else 0

    def test_no_floor_warning_with_the_default(self):
        self.assertEqual(self.floor_warnings(""), 0)

    def test_floor_warning_with_the_old_value(self):
        self.assertGreater(self.floor_warnings("material.min_amount = 1e-15"), 0)


if __name__ == "__main__":
    unittest.main()
