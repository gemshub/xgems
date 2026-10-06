# xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
#
# Checks that the solver modes (aia, sia, aop, sop, hop, shp) work on one example system, and that the
# warm start setting, the pH/Eh targets and the error handling behave as documented
# (see docs/sphinx/source/solver_guide.rst). The Optima modes are skipped when GEMS3K was
# built without Optima.

import os
import unittest

import numpy as np
from xgems import ChemicalEngine, ChemicalEngineDicts

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.join(HERE, "gems3k", "CASHNK1", "1TCi+CH-dat.lst")

OPTIMA = ChemicalEngine.builtWithOptima()
needs_optima = unittest.skipUnless(OPTIMA, "GEMS3K was built without Optima")

# GEMS3K status codes of a good result
OK = {"aia": 2, "sia": 6, "aop": 11, "sop": 15, "hop": 23, "shp": 27}


class TestSolverModes(unittest.TestCase):

    def setUp(self):
        self.engine = ChemicalEngine(PROJECT)
        self.T = self.engine.temperature()
        self.P = self.engine.pressure()
        # a copy: the array returned by elementAmounts() is a live view of the engine's amounts
        self.b = np.array(self.engine.elementAmounts())

    def solve(self, mode):
        self.engine.setSolverMode(mode)
        return self.engine.equilibrate(self.T, self.P, self.b.copy())

    def test_default_mode_is_aia(self):
        self.assertEqual(self.engine.solverMode(), "aia")
        self.assertEqual(self.engine.options.solver_mode, "aia")

    def test_sia_on_a_new_engine_starts_from_scratch(self):
        """With no previous result the first sia calculation runs as aia (status 2); the next one is warm."""
        self.assertEqual(self.solve("sia"), OK["aia"])
        self.assertEqual(self.solve("sia"), OK["sia"])

    def test_ipm_modes(self):
        for mode in ("aia", "sia"):
            with self.subTest(mode=mode):
                self.assertEqual(self.solve(mode), OK[mode])
                self.assertEqual(self.engine.solverMode(), mode)

    @needs_optima
    def test_optima_modes(self):
        for mode in ("aop", "sop"):
            with self.subTest(mode=mode):
                self.assertEqual(self.solve(mode), OK[mode])
                self.assertEqual(self.engine.solverMode(), mode)

    @needs_optima
    def test_hybrid_modes(self):
        self.solve("aia")  # previous result for the warm start
        for mode in ("hop", "shp"):
            with self.subTest(mode=mode):
                self.assertEqual(self.solve(mode), OK[mode])
                self.assertEqual(self.engine.solverMode(), mode)

    @needs_optima
    def test_modes_agree(self):
        """Every mode finds the same equilibrium of this system."""
        self.solve("aia")
        pH = self.engine.pH()
        G = self.engine.systemGibbsEnergy() if hasattr(self.engine, "systemGibbsEnergy") else None
        for mode in ("sia", "aop", "sop", "hop", "shp"):
            with self.subTest(mode=mode):
                self.solve(mode)
                self.assertAlmostEqual(self.engine.pH(), pH, delta=1e-3)
                if G is not None:
                    self.assertAlmostEqual(self.engine.systemGibbsEnergy(), G, delta=1e-3 * abs(G))

    def test_mode_name_is_case_insensitive(self):
        self.assertEqual(self.solve("AIA"), OK["aia"])
        self.assertEqual(self.engine.solverMode(), "aia")

    @needs_optima
    def test_warm_start_turns_cold_modes_warm(self):
        self.solve("aia")  # a warm mode needs a previous result to start from
        self.engine.setSolverMode("aia")
        self.engine.setWarmStart()
        self.assertEqual(self.engine.solverMode(), "sia")
        self.assertEqual(self.engine.equilibrate(self.T, self.P, self.b.copy()), OK["sia"])
        self.engine.setSolverMode("aop")
        self.assertEqual(self.engine.solverMode(), "sop")
        self.assertEqual(self.engine.equilibrate(self.T, self.P, self.b.copy()), OK["sop"])
        self.engine.setSolverMode("hop")
        self.assertEqual(self.engine.solverMode(), "shp")
        self.assertEqual(self.engine.equilibrate(self.T, self.P, self.b.copy()), OK["shp"])
        self.engine.setColdStart()
        self.assertEqual(self.engine.solverMode(), "hop")
        self.engine.setSolverMode("aop")
        self.assertEqual(self.engine.solverMode(), "aop")
        self.engine.setSolverMode("sop")
        self.engine.setColdStart()
        self.assertEqual(self.engine.solverMode(), "aop")

    def test_reequilibrate_with_warm_start(self):
        self.solve("aia")
        self.assertEqual(self.engine.reequilibrate(True), OK["sia"])
        self.assertEqual(self.engine.solverMode(), "sia")
        self.assertEqual(self.engine.reequilibrate(False), OK["aia"])
        self.assertEqual(self.engine.solverMode(), "aia")

    def test_set_through_options(self):
        opts = self.engine.options
        opts.solver_mode = "sia"
        self.engine.options = opts
        self.assertEqual(self.engine.solverMode(), "sia")

    def test_unknown_mode_is_rejected(self):
        for mode in ("rop", "native", "xx"):
            with self.subTest(mode=mode):
                with self.assertRaises(RuntimeError):
                    self.engine.setSolverMode(mode)
                self.assertEqual(self.engine.solverMode(), "aia")

    def test_rejected_options_keep_previous_mode(self):
        self.solve("aia")  # a warm mode needs a previous result to start from
        self.engine.setSolverMode("sia")
        opts = self.engine.options
        opts.solver_mode = "xx"
        with self.assertRaises(RuntimeError):
            self.engine.options = opts
        self.assertEqual(self.engine.solverMode(), "sia")
        self.assertEqual(self.engine.equilibrate(self.T, self.P, self.b.copy()), OK["sia"])

    @needs_optima
    def test_pH_target(self):
        self.engine.setSolverMode("aop")
        self.engine.setpHTarget(11.5)
        self.assertEqual(self.engine.equilibrate(self.T, self.P, self.b.copy()), OK["aop"])
        self.assertAlmostEqual(self.engine.pH(), 11.5, delta=1e-3)
        self.assertNotEqual(self.engine.controlConditionTitrant("pH"), 0.0)
        # the titrant stays in the engine's amounts, so start again from the original ones
        self.engine.clearControlConditions()
        self.engine.equilibrate(self.T, self.P, self.b.copy())
        self.assertGreater(self.engine.pH(), 12.0)
        self.assertEqual(self.engine.controlConditionTitrant("pH"), 0.0)

    @needs_optima
    def test_pH_and_Eh_targets_in_every_optima_mode(self):
        """aop, sop, hop and shp reach the same pH and Eh with the same titrant."""
        results = {}
        for mode in ("aop", "sop", "hop", "shp"):
            engine = ChemicalEngine(PROJECT)
            engine.setSolverMode("aia")
            engine.equilibrate(self.T, self.P, self.b.copy())  # previous result for the warm modes
            engine.setSolverMode(mode)
            engine.setpHTarget(11.5)
            engine.setEhTarget(0.2)
            self.assertEqual(engine.equilibrate(self.T, self.P, self.b.copy()), OK[mode])
            self.assertAlmostEqual(engine.pH(), 11.5, delta=1e-3)
            self.assertAlmostEqual(engine.Eh(), 0.2, delta=1e-3)
            results[mode] = (engine.controlConditionTitrant("pH"), engine.controlConditionTitrant("Eh"))
        for mode, (pH_titrant, Eh_titrant) in results.items():
            with self.subTest(mode=mode):
                self.assertAlmostEqual(pH_titrant, results["aop"][0], delta=1e-4)
                self.assertAlmostEqual(Eh_titrant, results["aop"][1], delta=1e-4)

    def test_ipm_modes_ignore_pH_target(self):
        if not OPTIMA:
            with self.assertRaises(RuntimeError):
                self.engine.setpHTarget(11.5)
            return
        self.engine.setpHTarget(11.5)
        self.solve("aia")
        self.assertGreater(self.engine.pH(), 12.0)
        self.assertEqual(self.engine.controlConditionTitrant("pH"), 0.0)

    def test_trace_regimes(self):
        self.solve("aia")
        regimes = self.engine.traceRegimes(of_interest=["Ca"])
        self.assertEqual(len(regimes), 1)
        self.assertEqual(regimes[0].name, "Ca")
        self.assertNotEqual(regimes[0].verdict, "FAILED")


class TestSolverModesDicts(unittest.TestCase):

    def setUp(self):
        self.engine = ChemicalEngineDicts(PROJECT)

    def test_modes(self):
        modes = ["aia", "sia"] + (["aop", "sop", "hop", "shp"] if ChemicalEngineDicts.builtWithOptima() else [])
        for mode in modes:
            with self.subTest(mode=mode):
                self.engine.setSolverMode(mode)
                self.assertEqual(self.engine.solverMode(), mode)
                self.assertTrue(self.engine.equilibrate().startswith("OK"))

    def test_built_with_optima_matches_engine(self):
        self.assertEqual(ChemicalEngineDicts.builtWithOptima(), OPTIMA)


if __name__ == "__main__":
    unittest.main(verbosity=2)
