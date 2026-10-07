# xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
#
# Checks that logging follows xgems.update_loggers(use_cout, logfile_name, log_level):
#  * the log files (ipmlog.txt, xGEMS.log, or the one named by the user) are created only when a message is written,
#  * the warnings of the Optima solver are logged at warning level, so they reach the console, a log file, or
#    nowhere, depending on update_loggers, instead of being printed to the console unconditionally.
#
# Each case runs in a fresh Python process in an empty directory, because the loggers are process-wide and the
# files are created in the current directory. The Optima cases need GEMS3K built with Optima and are skipped
# without it. The trigger is a pH target (25) that no equilibrium can reach, which makes the Optima line search
# warn that its Newton step is not a descent direction.

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

from xgems import ChemicalEngine

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.join(HERE, "gems3k", "CASHNK1", "1TCi+CH-dat.lst")

OPTIMA = ChemicalEngine.builtWithOptima()
needs_optima = unittest.skipUnless(OPTIMA, "GEMS3K was built without Optima")

PROLOGUE = f"""
import numpy as np
import xgems
from xgems import ChemicalEngine
engine = ChemicalEngine({PROJECT!r})
T, P, b = engine.temperature(), engine.pressure(), np.array(engine.elementAmounts())
"""

NORMAL_SOLVE = "engine.equilibrate(T, P, b.copy())\n"
UNREACHABLE_PH = """
engine.setSolverMode("aop")
engine.setpHTarget(25.0)
engine.equilibrate(T, P, b.copy())
"""


def run(code):
    """Run `code` in a new process in an empty directory.

    Returns (stdout, stderr, files created, {file: content}). The files are read after the process has ended,
    because the logger writes them out only then.
    """
    with tempfile.TemporaryDirectory() as cwd:
        proc = subprocess.run([sys.executable, "-c", textwrap.dedent(PROLOGUE + code)], cwd=cwd,
                              capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            raise AssertionError(f"subprocess failed:\n{proc.stderr}")
        files = sorted(os.listdir(cwd))
        contents = {f: open(os.path.join(cwd, f)).read() for f in files}
        return proc.stdout, proc.stderr, files, contents


class TestLogFiles(unittest.TestCase):

    def test_no_log_file_without_messages(self):
        # default logging, an ordinary solve that logs nothing: no file may appear
        out, err, files, contents = run(NORMAL_SOLVE)
        self.assertEqual(files, [])

    def test_no_log_file_when_logging_is_off(self):
        out, err, files, contents = run('xgems.update_loggers(False, "", 6)\n' + NORMAL_SOLVE)
        self.assertEqual(files, [])

    def test_user_log_file_is_created_on_first_message_only(self):
        out, err, files, contents = run('xgems.update_loggers(False, "mylog.txt", 3)\n' + NORMAL_SOLVE)
        self.assertNotIn("mylog.txt", files)


@needs_optima
class TestOptimaWarnings(unittest.TestCase):

    def test_warning_is_silent_when_logging_is_off(self):
        out, err, files, contents = run('xgems.update_loggers(False, "", 6)\n' + UNREACHABLE_PH)
        self.assertNotIn("OPTIMA WARNING", out + err)
        self.assertNotIn("Optima:", out + err)
        self.assertEqual(files, [])

    def test_warning_is_silent_above_warning_level(self):
        # level 4 = error: warnings are not logged
        out, err, files, contents = run('xgems.update_loggers(True, "log.txt", 4)\n' + UNREACHABLE_PH)
        self.assertNotIn("Optima:", out + err)
        self.assertNotIn("log.txt", files)

    def test_warning_reaches_the_console_when_asked(self):
        out, err, files, contents = run('xgems.update_loggers(True, "", 3)\n' + UNREACHABLE_PH)
        self.assertIn("Optima:", out + err)
        self.assertIn("not a descent direction", out + err)
        self.assertNotIn("***OPTIMA WARNING***", out + err)   # no longer printed directly by Optima

    def test_warning_goes_to_the_log_file_not_the_console(self):
        out, err, files, contents = run('xgems.update_loggers(False, "log.txt", 3)\n' + UNREACHABLE_PH)
        self.assertEqual((out + err).strip(), "")              # nothing on the console
        self.assertIn("log.txt", files)
        self.assertIn("Optima:", contents["log.txt"])
        self.assertIn("not a descent direction", contents["log.txt"])


if __name__ == "__main__":
    unittest.main()
