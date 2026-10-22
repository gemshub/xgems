Solver Guide
============

xGEMS has two solvers: the GEMS3K IPM solver and the **Optima** solver. This page shows how to choose a solver
and how to use it. Other pages give more detail:

* :doc:`advanced_usage`: batch runs, failed calls, absent elements and trace elements.
* :doc:`solver_settings`: the settings in the project files, and the known limits.
* :doc:`expert_settings`: settings for tests. Most users do not need them.

.. note::
   The Optima modes need a GEMS3K that is built with Optima. Check with ``ChemicalEngine.builtWithOptima()``.
   Without Optima, only the modes ``"aia"`` and ``"sia"`` work. Each other mode raises ``RuntimeError``.
   See :doc:`installation`.


The GEMS system files
---------------------

xGEMS reads a chemical system that GEM-Selektor exported. To export it, open a calculated equilibrium
(a SysEq record) in GEM-Selektor and calculate it again. Then choose *Data -> Write GEMS3K files*.
The `GEM-Selektor documentation <https://gemshub.github.io/site/start/gemselektor/documentation/>`_
explains how to make a project.

The export gives files that have the same name, for example ``MyTask``:

============================ =======================================================================================
File                         Content
============================ =======================================================================================
``MyTask-dat.lst``           The list of the other files. Give this file to xGEMS.
``MyTask-dch.json``          The chemical system: elements, species, phases and thermodynamic data.
``MyTask-ipm.json``          The solver settings (``pa_...``) and the parameters of the mixing models.
``MyTask-dbr-0-0000.json``   The bulk composition, the temperature and the pressure.
============================ =======================================================================================

All files must be in the same folder. The files are JSON (``-j``) or text (``-t``). The first entry of the
``-dat.lst`` file shows which one. The entries ``-f`` (JSON) and ``-o`` (text) show that a ThermoFun file
(``MyTask-fun.json``) holds the thermodynamic data.

.. code-block:: python

   from xgems import ChemicalEngine

   engine = ChemicalEngine("MyTask-dat.lst")
   b = engine.elementAmounts()
   print(engine.temperature(), engine.pressure())

Example systems are in the `tests/gems3k <https://github.com/gemshub/xgems/tree/master/tests/gems3k>`_ folder
of the xGEMS repository. The examples in this guide use ``CASHNK1/1TCi+CH-dat.lst``. This is a cement system
with an aqueous solution, a gas, a C-S-H solid solution and portlandite.


Select a solver
---------------

Each engine has one solver mode. The mode stays the same until you change it. The default is ``"aia"``.

.. code-block:: python

   engine.setSolverMode("aop")      # Optima solver, cold start
   engine.equilibrate(T, P, b)      # uses the selected mode
   print(engine.solverMode())       # shows the selected mode

======== ====================================================================================
Mode     What it does
======== ====================================================================================
aia      GEMS3K IPM solver, cold start. This is the default.
sia      GEMS3K IPM solver, warm start.
aop      Optima solver, cold start.
sop      Optima solver, warm start.
hop      GEMS3K IPM solver, cold start. Then Optima refines the result and checks it.
shp      The same as ``hop``, but with a warm start.
======== ====================================================================================

A cold start does not use a previous result. A warm start uses the previous result of the same engine.
A warm start costs less. A warm mode on a new engine runs as a cold start for the first call.

If Optima fails, the mode ``hop`` keeps the GEMS3K IPM result. Thus ``hop`` is not worse than ``aia``.

Optional: ``engine.setWarmStart()`` changes ``aia`` to ``sia``, ``aop`` to ``sop`` and ``hop`` to ``shp``.
``engine.setColdStart()`` reverses this. The options object has the same setting (``options.warmstart``).
You can also set the mode in the options object (``options.solver_mode``), in ``ChemicalEngineDicts`` and in C++:

.. code-block:: cpp

   xGEMS::ChemicalEngine engine("project-dat.lst");
   if (xGEMS::ChemicalEngine::builtWithOptima())
       engine.setSolverMode("aop");


Which mode to use
-----------------

.. list-table::
   :header-rows: 1
   :widths: 45 20 35

   * - Situation
     - Mode
     - Reason
   * - Most systems
     - ``aia``
     - It is the fastest solver on almost all systems.
   * - Solid solutions, melts, phase diagrams
     - ``aop``
     - Optima finds the correct phases more often.
   * - The ``aia`` mode fails or warns about the mass balance
     - ``aop``
     - Optima solves most of these cases.
   * - You want a check of the GEMS3K IPM result
     - ``hop``
     - Optima confirms or improves the result.
   * - Large systems (many species)
     - ``aia`` or ``hop``
     - Optima alone is much slower.
   * - Very little water
     - ``aia`` or ``hop``
     - Optima alone can fail when the water is almost gone.
   * - A sequence of similar points
     - ``sop`` after ``aop`` (or ``sia`` after ``aia``)
     - A warm start costs much less.
   * - A fixed pH or Eh
     - ``aop`` or ``sop``
     - The modes ``hop`` and ``shp`` also work.

A warm ``sia`` sequence can keep phases that are no longer stable after a phase boundary. Use ``aop`` and ``sop`` for such sequences.


Quick start
-----------

.. code-block:: python

   import numpy as np
   from xgems import ChemicalEngine

   engine = ChemicalEngine("CASHNK1/1TCi+CH-dat.lst")
   T, P = engine.temperature(), engine.pressure()
   b = np.array(engine.elementAmounts())         # a copy of the amounts

   engine.setSolverMode("aop")                   # cold start
   status = engine.equilibrate(T, P, b)
   print(engine.solverMode(), status, engine.pH(), engine.Eh())

   engine.setSolverMode("sop")                   # next, similar point: warm start
   status = engine.equilibrate(T + 5, P, b)

``equilibrate`` returns the GEMS3K status code. These codes show a good result:

======= ===== ===== ===== ===== ===== =====
Mode    aia   sia   aop   sop   hop   shp
======= ===== ===== ===== ===== ===== =====
Good    2     6     11    15    23    27
======= ===== ===== ===== ===== ===== =====

The next code in each group (3, 7, 12, 16, 24, 28) shows a result that you cannot fully trust.
The code after that (4, 8, 13, 17, 25, 29) shows a failure.
Two more codes are possible in each mode: 0 means that no calculation was necessary.
9 means a terminal error in GEMS3K. Restart the program after code 9.
``engine.converged()`` is true for a good result in all modes.

GEMS3K tries a failed call again in steps before it reports the failure. ``engine.rescueStep()`` shows
which step made the result. 0 means the call itself. See :ref:`batch-and-transport`.


See the results
---------------

``print(engine)`` writes the equilibrium state as tables. The tables show the temperature, the pressure,
the element amounts, the phases and the species. It works after each solver mode.

.. code-block:: python

   engine.equilibrate(T, P, b)
   print(engine)

``ChemicalEngine.equilibrate`` returns the status as a number. ``ChemicalEngineDicts.equilibrate`` returns the
status as text, for example ``OK after GEM calculation with LPP AIA``.


Fix the pH or Eh
----------------

The Optima modes (``aop``, ``sop``, ``hop`` and ``shp``) can fix the pH, the Eh (in volts), or both.
The solver adds the acid, base or electron donor that it needs. The targets stay active until you clear them.

.. code-block:: python

   engine.setSolverMode("aop")
   engine.setpHTarget(11.5)                      # optional second argument: tolerance
   engine.setEhTarget(0.2)                       # optional. You can use it alone.
   engine.equilibrate(T, P, b)

   print(engine.pH(), engine.Eh())
   print(engine.controlConditionTitrant("pH"))   # mol of titrant for the pH target

   engine.clearControlConditions()               # no targets

Note these points:

* The titrant becomes part of the bulk composition. The array from ``elementAmounts()`` is a live view,
  and it changes. Keep a copy (``b0 = np.array(engine.elementAmounts())``) and give the copy to ``equilibrate``.
* The modes ``aia`` and ``sia`` ignore the targets. They do not raise an error.
* A DBR file can give the targets in the optional fields ``pHtarget`` and ``Ehtarget``. GEMSGUI writes them
  when the system recipe has a target. The engine uses them only in the Optima modes. In other modes it ignores
  them and writes a warning. A value that is not a finite number is not used. A binary DBR file has no targets.
* A DBR target replaces a target that you set with ``setpHTarget()`` or ``setEhTarget()``. Your target comes back
  when a later DBR has no target. ``clearControlConditions()`` removes all targets.
* A target can be impossible, for example an Eh in a system with no redox couple. The status can be OK and the
  target can be missed. Compare ``engine.pH()`` and ``engine.Eh()`` with the target.
* A target very near a water stability line can be outside the water window of the model. The call then fails, correctly.


Advanced options
----------------

These options help in batch runs and transport codes. :doc:`advanced_usage` gives the details.

* **Warm start.** Solve the first point cold and the next points warm. A warm start costs much less.
  See :ref:`batch-and-transport`.
* **Rescue.** GEMS3K tries a failed call again in steps. ``engine.rescueStep()`` shows the step that made the result.
  Step 6 is a balanced result that is not converged. See :ref:`rescue-steps`.
* **Present phases.** Both solvers use the same rule to set an absent phase to 0. ``presentPhasesRel()`` gives the present phases.
* **Absent elements.** ``min_amount`` sets the placeholder amount of an element that is not in the recipe.
* **Trace elements.** ``traceRegimes`` shows if a trace element behaves linearly.

To change a solver setting, see :doc:`solver_settings`. The known limits of the solvers are on the same page.
