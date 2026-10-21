Advanced Usage
==============

This page is for batch runs, transport codes and special inputs. Read the :doc:`solver_guide` first.


.. _batch-and-transport:

Batch runs and reactive transport
---------------------------------

A batch run or a transport code solves many similar points. Each point needs a result. Use these rules:

* Solve the first point cold. Solve the next points warm (``sia``, ``sop`` or ``shp``).
* Give each call its own copy of the bulk composition. The array from ``elementAmounts()`` is a live view.
* Check each result. Do not stop the run on one failed point.

.. code-block:: python

   import numpy as np
   from xgems import ChemicalEngine

   engine = ChemicalEngine("project-dat.lst")
   T, P = engine.temperature(), engine.pressure()
   b0 = np.array(engine.elementAmounts())

   engine.setSolverMode("aop")                   # first point: cold start
   for i, x in enumerate(np.linspace(0.0, 1e-3, 50)):
       b = b0.copy()
       b[engine.indexElement("Ca")] += x         # the change of this point
       status = engine.equilibrate(T, P, b)
       if not engine.converged():
           print(i, "status", status, "rescue step", engine.rescueStep())
       engine.setSolverMode("sop")               # next points: warm start

**One engine for many nodes.** An engine holds one result. A warm start begins from the last result of that engine.
If one engine solves all nodes of a transport step, each node starts from the previous node.
This is good when neighbor nodes are similar. If they are very different, use a cold mode.
You can also use one engine for each group of similar nodes.


.. _rescue-steps:

Failed calls and the rescue
---------------------------

A call can end with an error or with a BAD result. Then GEMS3K tries the call again in steps.
It finishes each successful step at the exact composition, with the normal settings.
The setting is ``pa_FailRescue`` (default 2, on). ``engine.rescueStep()`` shows the step:

===== ====================================================================================
Step  What made the result
===== ====================================================================================
0     The call itself. No rescue was necessary.
1     A warm restart from the previous result.
2     A composition that is changed by a very small amount.
3     A walk in equal steps from the previous result to the new composition.
4     Looser settings.
5     The other solver.
6     No step succeeded. You get the best attempt that keeps the element totals. Its status is always BAD.
===== ====================================================================================

Keep ``pa_FailRescue = 2`` for transport codes. A transport code needs a result for each node.
The rescue costs more iterations only on calls that fail.

Two limits keep the cost of a call that no step can save small:

* ``pa_RescueBudget`` (default 100000). After the call has used this number of iterations, the rescue skips
  steps 1 to 4. It still tries the other solver (step 5) one time.
* ``pa_OptimaStallNetFactor`` (default 5). When an Optima solve stops making progress, Optima solves one more time
  with the progress check off. This last solve gets at most 5 times the iterations of the solve that stopped.

A result from step 6 keeps the mass balance, but it is not converged. Mark this point, or solve it again with a smaller change.
Use ``pa_FailRescue = 1`` if you prefer an error to a BAD result. Use ``pa_FailRescue = 0`` to switch the rescue off.

.. code-block:: python

   status = engine.equilibrate(T, P, b)
   step = engine.rescueStep()
   if step == 6:
       bad_points.append(i)                      # balanced, but not converged
   elif step > 0:
       rescued_points.append((i, step))          # a good result after a rescue

``ChemicalEngineDicts`` gives the step as a property (``engine.rescue_step``). In C++, use ``engine.rescueStep()``.

.. _suggested-settings:

Suggested settings for common cases
-----------------------------------

Most projects need no change. The table gives the settings that help in some common cases.
They are suggestions for one project or one run. They do not change the defaults.
The GEMS3K `advanced usage <https://github.com/gemshub/GEMS3K/blob/develop_optima2/Docs/SOLVER-ADVANCED-USAGE.md>`_ page has the same table.
Put them in the ``-ipm`` file of the project (see :ref:`project-file-settings`). The code after the table changes them from Python.

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Case
     - Settings and mode
     - Why
   * - Reactive transport and batch runs: each node must have a result
     - Keep the defaults ``pa_FailRescue = 2`` and ``pa_RescueBudget = 100000``. Solve the first point cold and the next
       points warm (``sia``, ``sop``). ``pa_PSM = 0`` for large runs.
     - The rescue tries a failed node again in steps. The limit stops a node that no step can save.
       ``pa_PSM = 0`` writes no log file.
   * - Some calls run for a very long time
     - ``pa_OptimaCallCap = -1``
     - Sets a cap from the size of the system. The results do not depend on the computer.
       ``pa_OptimaMaxSeconds`` also limits each solve, but then the results depend on the speed of the computer.
   * - Phase diagrams and Pourbaix maps where a failed point is acceptable
     - ``pa_FailRescue = 0``. Use ``aop`` or ``sop``, not ``sia``, along a row.
     - A point that cannot be solved, for example a target outside the water window of the model, then costs much
       less. You get the error instead of the best attempt. ``sia`` can keep phases after a phase boundary.
   * - Water used up by design (for example cement hydration with little water)
     - ``pa_OptimaDryOut = 1``. If ``aop`` or ``sop`` still fails: ``pa_OptimaAcceptRepair = 1`` and ``pa_GAS = 0.002``,
       or use ``aia`` or ``hop``.
     - Optima solves the solids without the water phase when the water runs out. The pH of a result without water
       is reported as 0.
   * - No redox couple in the system
     - ``pa_OptimaRedoxFree = 2``, or give an Eh target.
     - The bulk composition does not fix the Eh. The setting aims at pe 4.
   * - Solid solutions with end-members that are almost the same
     - ``pa_DK = 1e-7``
     - Gives a clearly more accurate composition with the GEMS3K IPM solver.
   * - Sorption models
     - ``pa_DGC = 0.001``
     - Damps the non-ideal corrections between iterations.
   * - Absent phases must show exactly zero
     - ``pa_OptimaZeroAbsent = 1``
     - Optima keeps an absent species at a very small floor amount by default.
   * - Large systems (hundreds of species)
     - Mode ``aia``, or ``hop`` for a checked result. With Optima: ``pa_OptimaTol`` a little larger.
     - Optima alone is much slower on large systems.
   * - Trace elements
     - A stricter ``pa_DHB``. Put the element in the recipe, or set ``min_amount``.
     - Keeps the mass balance of small amounts exact.
   * - Find the cause of a failed calculation
     - ``pa_PSM = 2`` or ``3``
     - Writes warnings or a detailed trace to ``ipmlog.txt``.

xGEMS has no call that changes a setting after the project is loaded. Copy the project and edit the copy.
This code works for JSON projects:

.. code-block:: python

   import json, pathlib, shutil, tempfile
   from xgems import ChemicalEngine

   def project_with_settings(lst, **settings):
       """Copy a JSON project to a new folder and set pa_ values in its -ipm file."""
       src, dst = pathlib.Path(lst).parent, pathlib.Path(tempfile.mkdtemp())   # the folder is outside the project
       shutil.copytree(src, dst, dirs_exist_ok=True)
       for ipm in dst.glob("*-ipm.json"):
           data = json.loads(ipm.read_text())
           data[0]["ipm"].update(settings)       # replaces a value, or adds a setting
           ipm.write_text(json.dumps(data, indent=2))
       return str(dst / pathlib.Path(lst).name)

   lst = project_with_settings("project-dat.lst", pa_PSM=0, pa_OptimaCallCap=-1)
   engine = ChemicalEngine(lst)


Present phases
--------------

Optima never sets an amount to exactly zero. A phase that does not form stays at a very small floor amount.
A fixed limit in mol can count this phase as present.

``presentPhasesRel()`` uses the GEMS3K rule. A phase is present when it holds more than 1e-9 of the total amount of all phases.
With this rule, the Optima modes and the GEMS3K IPM modes show the same phases.

.. code-block:: python

   print(engine.presentPhasesRel())              # ChemicalEngine
   print(engine.present_phases_rel(1e-9))        # ChemicalEngineDicts


Elements that are not in the recipe
-----------------------------------

Each element of the chemical system needs an amount larger than zero. If the bulk composition does not contain an element,
xGEMS gives it a small placeholder amount, ``min_amount``. The default is **1e-11 mol**.
The charge element ``Zz`` is always 0.

``min_amount`` is a property of ``Material``. It is also an argument of these methods of ``ChemicalEngineDicts``:
``equilibrate``, ``setB``, ``clear``, ``set_bulk_composition`` and ``get_b_from_formula``.
In C++ the default is ``xGEMS::default_min_element_amount``.

.. code-block:: python

   mat = Material(engine, "recipe")
   mat.add("H2O", 1.0, "kg")
   mat.add({"Ca": 1e-3, "Si": 1e-3, "O": 3e-3})     # N and other elements get 1e-11 mol
   mat.min_amount = 1e-9                            # another placeholder amount
   mat.min_amount = 0.0                             # no placeholder

**Why not a smaller value.** The Optima modes keep each species above a small floor amount.
If the total amount of an element is smaller than its species need at this floor, the mass balance cannot be correct.
GEMS3K then repairs the mass balance after the solve. This costs time, and GEMS3K writes two warnings:

* ``Optima: N element(s) have less material than the solver's floor amount allows ...``. It names the elements and the amount that is enough.
* ``Mass balance is still off after repair (worst element X, ...x its tolerance)``.

The modes ``aia`` and ``sia`` do not use the floor. They do not give the first warning.
The `GEMS3K logging documentation <https://github.com/gemshub/GEMS3K/blob/master/Docs/spdlog-doc.md>`_ explains where the warnings go.

**If the warnings continue.** Some systems need a larger ``min_amount``. The first warning gives the amount to use.
Raise ``min_amount`` only as far as necessary. The placeholder is real material.
For example, 1e-7 mol of an element in 1 kg of water is a concentration of 0.1 micromolal.

**For trace chemistry.** If 1e-11 mol of an extra element is too much (for example a radionuclide at ppb level),
put the element in the recipe. You can also set ``min_amount`` to 1e-15 or 0.


Trace elements
--------------

``traceRegimes`` shows if a trace element behaves linearly. A linear element has a distribution that scales with its amount.
The phases do not change.

The method solves the system cold again in the current solver mode. It scales the trace amounts by each factor
and compares the phase fractions. The engine state is restored afterward.

.. code-block:: python

   engine.equilibrate(T, P, b)
   for r in engine.traceRegimes(of_interest=["Ca", "Si"]):
       print(r.name, r.verdict, r.phases)

An element is a trace element when its amount is at most ``trace_rel`` (default 1e-6) times the total.
Use ``of_interest`` to check other elements. ``tol`` (default 1e-3) is the largest change of a phase fraction that is still linear.

========= ==================================================================================
Verdict   Meaning
========= ==================================================================================
LINEAR    The phase fractions do not change at any factor.
SATURATED A pure phase of the element is present. The dissolved amount is constant.
BOUNDARY  A pure phase of the element appears or disappears (see ``boundaryPhases``).
NONLINEAR The distribution changes with the amount. No phase appears or disappears.
STRANDED  All of the element is in one mixed phase that has only a trace amount.
          The verdict can have two parts, for example ``STRANDED+NONLINEAR``.
FAILED    One of the new solves failed.
========= ==================================================================================

``ChemicalEngineDicts`` returns ``{element name: TraceRegime}``.
