Solver Guide
============

This page describes the two equilibrium solvers xGEMS can use: the GEMS3K IPM solver, which has always been
there, and the **Optima** solver, which is new in GEMS3K 5.0 and xGEMS 2.2. It explains
when to use which, how to switch from Python and C++, how the classic IPM solver works and how to tune it,
how to fix the pH or Eh with Optima, how to check trace elements, and what to expect from each solver.

.. note::
   The Optima modes need a GEMS3K built with Optima. Check with
   ``ChemicalEngine.builtWithOptima()``. Without Optima only the ``"aia"`` and ``"sia"`` modes are available, and
   every other call below raises ``RuntimeError``. See :doc:`installation`.


The GEMS system files
---------------------

xGEMS works on a chemical system exported from GEM-Selektor (or created by the GEMS3K code). How to create a
modeling project in GEM-Selektor is described in the `GEM-Selektor documentation
<https://gemshub.github.io/site/start/gemselektor/documentation/>`_. To export it, open a calculated equilibrium
(a SysEq record), re-calculate it, and choose *Data -> Write GEMS3K files*. The export gives a small
set of files that share a name, for example ``MyTask``:

============================ =======================================================================================
File                         What it holds
============================ =======================================================================================
``MyTask-dat.lst``           The list of the three files below. This is the file you give to xGEMS.
``MyTask-dch.json``          The chemical system: elements, species and phases, and their thermodynamic data.
``MyTask-ipm.json``          The solver settings (the ``pa_...`` entries) and the parameters of the mixing models.
``MyTask-dbr-0-0000.json``   The bulk composition, temperature and pressure of a calculation. When it was written
                             after a calculation, it also holds the equilibrium result.
============================ =======================================================================================

The files come as key-value text (``.dat``) or as JSON (``.json``). The first entry of ``-dat.lst`` tells which one:
``-t`` for text and ``-j`` for JSON. The three files must be in the same folder as the ``-dat.lst`` file.

There is also a ThermoFun option. With ``-f`` (JSON) or ``-o`` (text) as the first entry, the ``-dat.lst`` lists one
more file, ``MyTask-fun.json``. The thermodynamic data of the substances then sits in that ThermoFun file, and
the thermodynamic properties are calculated with ThermoFun at the temperature and pressure
of each calculation, instead of being read from the tables in the ``-dch`` file. The kaolinite examples below use it.

The solver settings live in the ``-ipm`` file; they are explained in :ref:`project-file-settings`.

.. code-block:: python

   from xgems import ChemicalEngine

   engine = ChemicalEngine("MyTask-dat.lst")      # reads the three files
   b = engine.elementAmounts()                    # the bulk composition from the -dbr file
   print(engine.temperature(), engine.pressure())

Ready-made example systems are in the `tests/gems3k <https://github.com/gemshub/xgems/tree/master/tests/gems3k>`_ folder of the xGEMS repository.
The examples in this guide use ``CASHNK1/1TCi+CH-dat.lst``, a cement system (Ca, Si, N, H, O) with an aqueous
solution, a gas, a C-S-H solid solution and portlandite. The folders are:

* `CASHNK1 <https://github.com/gemshub/xgems/tree/master/tests/gems3k/CASHNK1>`_ and `CASHNK1-json <https://github.com/gemshub/xgems/tree/master/tests/gems3k/CASHNK1-json>`_: the system as JSON files (the two
  folders hold the same files)
* `CASHNK1-keyvalue <https://github.com/gemshub/xgems/tree/master/tests/gems3k/CASHNK1-keyvalue>`_: the same system as text files
* `CASHNK1-T-json <https://github.com/gemshub/xgems/tree/master/tests/gems3k/CASHNK1-T-json>`_: JSON files with a finer temperature grid (1 K steps from 274.15 K,
  instead of 5 K steps)
* `Kaolinite-funjson <https://github.com/gemshub/xgems/tree/master/tests/gems3k/Kaolinite-funjson>`_, `Kaolinite-funkeyvalue <https://github.com/gemshub/xgems/tree/master/tests/gems3k/Kaolinite-funkeyvalue>`_ and
  `Kaolinite-T-funjson <https://github.com/gemshub/xgems/tree/master/tests/gems3k/Kaolinite-T-funjson>`_: a second system, a kaolinite pH titration (``pHtitr``),
  in the same formats

Download a folder, and give the path of its ``-dat.lst`` file to ``ChemicalEngine``.

Selecting a solver
------------------

The solver is chosen per engine, with one call, and stays selected until you change it. The default is
``"aia"``, the GEMS3K IPM solver from scratch, so nothing needs to be set to get it:

.. code-block:: python

   from xgems import ChemicalEngine

   engine = ChemicalEngine("project-dat.lst")

   engine.setSolverMode("aia")      # GEMS3K IPM solver, start from scratch (default)
   engine.setSolverMode("sia")      # GEMS3K IPM solver, start from the previous result
   engine.setSolverMode("aop")      # Optima solver, start from scratch
   engine.setSolverMode("sop")      # Optima solver, start from the previous result
   engine.setSolverMode("hop")      # IPM solver, then Optima checks the answer
   engine.setSolverMode("shp")      # the same, starting from the previous result

   engine.equilibrate(T, P, b)      # uses the selected solver
   print(engine.solverMode())       # which one is selected

The same from the options object, from the dictionary interface and from C++. The Optima modes (``"aop"``, ``"sop"``,
``"hop"`` and ``"shp"``) need a GEMS3K built with Optima: without it they raise ``RuntimeError`` in every one of these forms, while ``"aia"``
and ``"sia"`` always work. Check ``ChemicalEngine.builtWithOptima()`` when your code must run on either:

.. code-block:: python

   if ChemicalEngine.builtWithOptima():     # "aop" and "sop" need a GEMS3K built with Optima
       opts = engine.options
       opts.solver_mode = "aop"
       engine.options = opts                # RuntimeError for an unknown mode; the engine keeps its previous mode

   from xgems import ChemicalEngineDicts
   engine = ChemicalEngineDicts("project-dat.lst")
   engine.setSolverMode("aia")              # always available
   if ChemicalEngineDicts.builtWithOptima():
       engine.setSolverMode("aop")

.. code-block:: cpp

   xGEMS::ChemicalEngine engine("project-dat.lst");
   if (xGEMS::ChemicalEngine::builtWithOptima())
       engine.setSolverMode("aop");


The two solvers
---------------

**GEMS3K IPM solver (modes** ``"aia"`` **and** ``"sia"`` **).** An interior-point method followed by a
refinement of the mass balance. It is fast on almost every system and is the default (``"aia"``).

**Optima solver (modes** ``"aop"`` **and** ``"sop"`` **).** A second solver, based on the open-source
`Optima <https://github.com/gemshub/optima>`_ optimisation library. It is more robust for solid solutions,
melts and other systems where several mixing phases compete, and it is the solver that fixes
the pH or Eh.

**Hybrid (modes** ``"hop"`` **and** ``"shp"`` **).** The GEMS3K IPM solver runs first, then Optima refines
and checks its answer. If the Optima part fails, the IPM answer is kept and the result is flagged, so a hybrid
calculation is never worse than the IPM solver alone.

======== ============ ==========================================================================
Mode     GEMS3K mode  What it does
======== ============ ==========================================================================
aia      AIA          GEMS3K IPM solver. Builds its own starting point from scratch (cold start).
sia      SIA          GEMS3K IPM solver. Starts from the previous result (warm start).
aop      AOP          Optima solver, starting from scratch (cold start).
sop      SOP          Optima solver, starting from the previous result (warm start).
hop      HOP          IPM solver from scratch, then Optima refines and checks the answer.
shp      SHP          As ``hop``, but the IPM part starts from the previous result.
======== ============ ==========================================================================

**Cold and warm start.** The mode name already says how the calculation starts. The warm start
setting does the same from the other side: when it is on (``engine.setWarmStart()``,
``options.warmstart = True`` or ``engine.reequilibrate(True)``), ``"aia"`` runs as ``"sia"``, ``"aop"`` as
``"sop"`` and ``"hop"`` as ``"shp"``. ``engine.setColdStart()`` turns every warm mode cold again, and ``engine.solverMode()`` returns the
mode that actually runs.

A warm mode needs a previous result to start from. On a new engine the first ``"sia"`` calculation therefore runs
from scratch and returns the AIA status code (2); the following ones are warm (status 6).



Which mode to use
-----------------

.. list-table::
   :header-rows: 1
   :widths: 40 18 42

   * - Situation
     - Mode
     - Why
   * - Most systems, single calculations
     - ``aia``
     - Fastest on almost every system, often by a large margin.
   * - Many species (hundreds or more)
     - ``aia``, or ``hop``
     - Optima alone is much slower on large systems; ``hop`` adds Optima's check for little more than the IPM cost.
   * - Solid solutions, melts, miscibility gaps, phase diagrams
     - ``aop``
     - Optima finds the correct set of stable phases far more often where several mixing phases compete.
   * - The GEMS3K IPM solver fails or warns about the mass balance
     - ``aop``
     - Optima solves most of these cases.
   * - You want the IPM answer checked
     - ``hop``
     - Optima confirms or improves it, and the IPM answer is kept if Optima fails.
   * - Sweeps in T or composition, titrations, transport steps
     - ``sop`` after ``aop`` (or ``sia`` after ``aia``, ``shp`` after ``hop``)
     - Starting from the previous result is many times cheaper.
   * - Sequences that cross phase boundaries
     - ``aop`` or ``sop``
     - Warm-started ``sia`` can keep phases that stopped being stable.
   * - Fixing the pH or Eh
     - ``aop`` or ``sop``
     - ``aop`` or ``sop`` (``hop`` and ``shp`` also support it).
   * - Systems with very little water
     - ``aia`` or ``hop``
     - Optima alone can fail when water is close to running out.


Quick start (Python)
--------------------

.. code-block:: python

   import numpy as np
   from xgems import ChemicalEngine

   engine = ChemicalEngine("CASHNK1/1TCi+CH-dat.lst")
   print(ChemicalEngine.builtWithOptima())       # True when Optima is available

   T, P = engine.temperature(), engine.pressure()
   b = np.array(engine.elementAmounts())         # a copy, see the note on pH/Eh control

   engine.setSolverMode("aop")                   # Optima, cold start
   status = engine.equilibrate(T, P, b)
   print(engine.solverMode(), status, engine.pH(), engine.Eh())

   engine.setSolverMode("sop")                   # next, similar point: warm start
   status = engine.equilibrate(T + 5, P, b)

The return value of ``equilibrate`` is the GEMS3K status code. The good ones are 2 (AIA), 6 (SIA), 11 (AOP), 15 (SOP),
23 (HOP) and 27 (SHP). The next code in each group (3, 7, 12, 16, 24 and 28) means the result is not fully
trustworthy, and the one after that (4, 8, 13, 17, 25 and 29) means failure.


Looking at the results
----------------------

``print(engine)`` writes the whole equilibrium state as tables: temperature and pressure, the element amounts,
the composition and properties of every phase, and the amount, activity and chemical potential of every species.
It works on a ``ChemicalEngine`` and on a ``ChemicalEngineDicts``, after any solver mode:

.. code-block:: python

   engine.equilibrate(T, P, b)
   print(engine)

The first tables of the output for the example system look like this (the species table, with all species,
follows):

.. code-block:: text

   ==============================================================================
   Temperature[K]           Temperature[C]           Pressure[MPa]
   ------------------------------------------------------------------------------
   298.15                   25                       0.1
   ==============================================================================
   Element                  InputAmount[mol]
   ------------------------------------------------------------------------------
   Ca                       4.000000e-01
   H                        1.118167e+02
   Nit                      2.000000e+00
   O                        5.671037e+01
   Si                       2.000000e-01
   Zz                       0.000000e+00
   ==============================================================================
   Property                 aq_gen                   gas_gen                  CSHK
   ------------------------------------------------------------------------------
   PhaseAmount[mol]         5.551774e+01             1.000350e+00             9.027775e-02
   PhaseMass[kg]            1.000666e+00             2.802720e-02             3.356960e-02
   PhaseVolume[m^3]         1.001646e-03             2.479840e-02             1.291614e-05
   ...

The return value of ``equilibrate`` is separate from this. ``ChemicalEngine.equilibrate`` returns the status as a
number, so ``print(engine.equilibrate(T, P, b))`` prints ``2``. ``ChemicalEngineDicts.equilibrate`` returns it as
text, so ``print(engine.equilibrate())`` prints ``OK after GEM calculation with LPP AIA``. The text names the solver
that ran, for example ``OK after GEM calculation via Optima with cold initial approximation (AOP)``.


Fixing the pH or Eh
-------------------

With an Optima mode (``aop``, ``sop``, ``hop`` or ``shp``) you can ask for a given pH, a given Eh (in volts), or both. The solver adds whatever
acid, base or electron donor the target needs, and reports how much it added. The targets stay active
for every later calculation until you clear them.

.. code-block:: python

   engine.setSolverMode("aop")
   engine.setpHTarget(11.5)                      # optional second argument: tolerance
   engine.setEhTarget(0.2)                       # optional; may be used alone
   engine.equilibrate(T, P, b)

   print(engine.pH(), engine.Eh())
   print(engine.controlConditionTitrant("pH"))   # mol of titrant solved for the pH target
   print(engine.controlConditionTitrant("Eh"))

   engine.clearControlConditions()               # back to an unconstrained calculation

Things to know:

* **The titrant becomes part of the bulk composition.** It is added to the engine's own amounts, and the
  array returned by ``elementAmounts()`` is a live view of them, so it changes too. To start the next point
  from the original composition, keep a copy (``b0 = np.array(engine.elementAmounts())``) and pass ``b0`` to
  ``equilibrate``. Passing the view again would keep the pH at the target after you clear the conditions.
  The dictionary interface is not affected, because you pass a new dictionary each time.
* **The ``aia`` and ``sia`` modes ignore the targets without an error.** The pH is not constrained and
  ``controlConditionTitrant`` returns 0. Select an Optima mode first.
* A target may not be reachable, for example Eh in a system with no redox couple, or a pH far outside the
  physical range. The status can then still say OK while the target was missed (seen with ``sop``, ``hop`` and ``shp``
  at pH 25). After the call, compare ``engine.pH()`` and ``engine.Eh()`` with the target, and do not rely on the Eh
  of a system without a redox couple.


Trace elements
--------------

``traceRegimes`` tells you whether a trace element behaves linearly: that is, whether scaling its amount
only scales its distribution, or changes which phases are present. It re-solves the system cold,
in the current solver mode, with the trace amounts scaled by each factor, and compares the phase fractions.
The engine state is restored afterwards.

.. code-block:: python

   engine.equilibrate(T, P, b)
   for r in engine.traceRegimes(of_interest=["Ca", "Si"]):
       print(r.name, r.verdict, r.phases)

An element is trace when its amount is at most ``trace_rel`` (default 1e-6) times the total; naming
elements in ``of_interest`` checks them whatever their amount. ``tol`` (default 1e-3) is the largest
phase-fraction change still counted as linear. The verdicts are:

========= ==================================================================================
LINEAR    Phase fractions are unchanged at every factor.
SATURATED A pure phase of the element is present and its dissolved amount is constant.
BOUNDARY  A pure phase of the element appears or disappears (see ``boundaryPhases``).
NONLINEAR The distribution changes with the amount, with no phase appearing or disappearing.
STRANDED  All of the element sits in one multi-species phase present only in a trace amount.
          May be combined, for example ``STRANDED+NONLINEAR``.
FAILED    One of the re-solves failed.
========= ==================================================================================

The dictionary interface returns ``{element name: TraceRegime}``.


GEMS3K IPM solver (aia and sia modes)
-------------------------------------

The IPM solver is the classic GEMS3K algorithm and the default of xGEMS. It minimises the Gibbs energy of the
system by an interior-point method (IPM) and has these stages:

1. **Starting point.** With AIA (automatic initial approximation) it builds its own starting point from
   scratch by a linear-programming step. With SIA (smart initial approximation, ``warmstart`` on) it
   starts from the previous result held in the engine, which is much cheaper for small changes.
2. **Main iteration.** Interior-point iterations on the species amounts, recalculating the activity
   coefficients of the non-ideal phases as it goes (``pa_PD``), until the convergence criterion ``pa_DK`` is met or
   ``pa_IIM`` iterations are used.
3. **Mass-balance refinement.** A separate step makes the amount of every element in the answer match the
   amount put in, within ``pa_DHB`` (and ``pa_DP`` iterations).
4. **Phase selection.** Phases that should be present but were left out are added back, and phases that are
   no longer stable are removed (``pa_PC``, ``pa_DF``, ``pa_DFM``), and the calculation repeats.
5. **Clean-up.** Tiny leftover amounts are removed (``pa_PRD``, ``pa_GAS``).

When AIA does not converge, GEMS3K retries automatically, from a slightly changed composition
(``pa_ColdRetryNudges``) and then from scratch, before it reports a failure.

The IPM solver is the fastest on almost every system. Its weak point is strongly non-ideal systems where
several solid solutions or melts compete: it can then leave out a phase that should be present, or
keep one that should not. Use ``"aop"`` there.

The GEMS3K IPM solver writes diagnostics to ``ipmlog.txt`` in the working directory; ``pa_PSM`` controls how
much (0 = nothing, 2 = also warnings, 3 = detailed trace). The log level of the xGEMS and GEMS3K messages on
screen is set with ``update_loggers``; what each of its values means, and which loggers exist, is described in the
`GEMS3K logging documentation <https://github.com/gemshub/GEMS3K/blob/master/Docs/spdlog-doc.md>`_.


.. _project-file-settings:

Settings in project files
-------------------------

Both solvers read their settings from the project's ``-ipm`` file (``pa_...`` entries), once, when the
project is loaded. xGEMS has no call to change them afterwards: edit the file and create a new engine.
Settings missing from the file keep their default. Files exported by earlier versions do not contain the
newer settings, so add them yourself when you need another value.

To add a setting to a JSON project, put it inside the ``"ipm"`` block of ``…-ipm.json``:

.. code-block:: json

   "pa_PLLG": 10000,
   "pa_OptimaTol": 1e-09,
   "pa_OptimaMaxSeconds": 60,

In a text project (``…-ipm.dat``), add it after the ``<END_DIM>`` line, for example
``<pa_OptimaTol>  1e-09``. Keep each setting only once in the file: if a name appears twice, the later value
is used. A setting placed before ``<END_DIM>`` in a text file is not read.

Some settings, such as ``pa_PE``, ``pa_DHB`` and ``pa_DG``, apply to both solvers.


GEMS3K IPM solver settings
~~~~~~~~~~~~~~~~~~~~~~~~~~

These settings are in most project files already.

.. list-table::
   :header-rows: 1
   :widths: 17 8 45 30

   * - Setting
     - Default
     - What it does and when it is useful
     - Other values
   * - ``pa_DK``
     - 1e-6
     - How precisely the GEMS3K IPM solver must converge. Smaller gives a more precise answer at the cost of more iterations. Solid solutions whose end-members are nearly alike: 1e-7 gives a noticeably more accurate composition.
     - Any small positive number.
   * - ``pa_DHB``
     - 1e-13
     - How closely the amount of each element in the answer must match the amount put in, as a fraction of that amount. Loosen it slightly if a calculation fails only because of this check; tighten it for systems with trace elements.
     - Usually between 1e-9 (looser) and 1e-15 (stricter).
   * - ``pa_DT``
     - 0
     - How the check above is applied. 0 applies it the same way to every element. Systems mixing major and trace elements where the relative check is too strict for the major ones.
     - A value of -6 or less = for major elements, use a fixed amount instead of a fraction.
   * - ``pa_IIM``
     - 7000
     - Maximum number of iterations of the main calculation before it gives up. Raise it for large or difficult systems that run out of iterations.
     - Up to 9999.
   * - ``pa_DP``
     - 130
     - Maximum number of iterations of the step that makes the element totals add up. Raise it if that step runs out of iterations on a difficult system.
     - Any positive number.
   * - ``pa_DW``
     - 1 (on)
     - Treats running out of the iterations above as an error. Turn off only to inspect a result that would otherwise be rejected.
     - 0 = do not treat it as an error.
   * - ``pa_PE``
     - 1 (on)
     - Requires the answer to be electrically neutral. Leave on for any system with ions; off only for systems with no charged species.
     - 0 = off.
   * - ``pa_PC``
     - 2
     - How the solver decides which phases are present. 2 is the current method, which also adds back phases that were wrongly left out. Leave at 2. The old method is kept for comparison with old results.
     - 1 = the old method.
   * - ``pa_DF``
     - 0.01
     - How clearly a missing phase must be stable before it is added to the answer. Lower it if phases that should appear are missed; raise it if phases flicker in and out along a sweep.
     - Any small positive number; smaller adds phases more readily.
   * - ``pa_DFM``
     - 0.01
     - How clearly a present phase must be unstable before it is removed from the answer. Adjust together with ``pa_DF`` when phases flicker in and out along a sweep.
     - Any small positive number.
   * - ``pa_PD``
     - 2
     - How often activity coefficients (the corrections for non-ideal behaviour) are recalculated. 2 = at every iteration. Leave at 2. Lower values can save time on nearly ideal systems.
     - 0 = once at the start; 1 = only while fixing the element totals; 3 = only in the main calculation.
   * - ``pa_AG``
     - 1
     - Damping of the non-ideal corrections between iterations. Lower it if a strongly non-ideal system oscillates and fails to converge.
     - From -1 to 1.
   * - ``pa_DGC``
     - 0
     - A second damping setting that works together with the one above. Use 0.001 for systems with sorption.
     - From -1 to 1.
   * - ``pa_PRD``
     - -5
     - A final clean-up that removes tiny amounts of species left over at the end of the calculation. -5 means amounts below 1e-5 are cleaned up. Make it more negative if real trace amounts are being removed; 0 to see the raw result.
     - 0 = off; -6 or less = clean up only smaller amounts.
   * - ``pa_GAS``
     - 0.001
     - How strict the final clean-up is. A larger value keeps more trace amounts in the answer. Systems with very little water (0.002 together with ``pa_OptimaAcceptRepair``), or when trace phases matter.
     - Any small positive number.
   * - ``pa_DS``
     - 1e-20
     - Smallest amount of a phase, in moles, that is still reported as present. Raise it to hide meaningless trace phases in the output.
     - Any small positive amount.
   * - ``pa_XwMin``
     - 1e-13
     - If the amount of water falls below this (moles), the aqueous solution is removed from the answer. Drying or nearly dry systems, to control when the aqueous solution disappears.
     - Any small positive amount.
   * - ``pa_ScMin``
     - 1e-13
     - If the amount of a solid that carries a sorption surface falls below this (moles), the sorption phase is removed. Sorption systems where the sorbent dissolves almost completely.
     - Any small positive amount.
   * - ``pa_PhMin``
     - 1e-20
     - If a solution phase other than the aqueous one falls below this amount (moles), it is removed with all its species. Systems with many solid solutions or melts present in trace amounts.
     - Any small positive amount.
   * - ``pa_DcMin``
     - 1e-33
     - If a species inside a mixed phase falls below this amount (moles), it is removed. Rarely needs changing.
     - Any small positive amount.
   * - ``pa_DB``
     - 1e-17
     - Smallest amount of an element, in moles, that the solver works with in the bulk composition (charge excluded). Systems with elements at extremely low amounts.
     - Any small positive amount.
   * - ``pa_ICmin``
     - 1e-5
     - Below this ionic strength (molal), activity coefficients of aqueous species are taken as 1, i.e. the solution is treated as ideal. Very dilute waters, if the ideal treatment starts too early or too late.
     - Any small positive number.
   * - ``pa_DG``
     - 1000
     - The solver internally scales the whole system to this total number of moles. Results are reported in the original amounts. Very small or very large systems; switch off only to compare with unscaled results.
     - Below 1e-4 = no scaling.
   * - ``pa_EPS``
     - 1e-10
     - How precisely the automatic starting point (used by AIA) is computed. Loosen it if the starting point cannot be found; tighten it for systems with trace elements.
     - From 1e-6 (looser) to 1e-14 (stricter).
   * - ``pa_DFYw``, ``pa_DFYaq``, ``pa_DFYid``, ``pa_DFYr``, ``pa_DFYh``, ``pa_DFYc``
     - 1e-5
     - Small amounts (moles) given to water, aqueous species, species of ideal and non-ideal mixed phases, and pure phases that are zero in the automatic starting point, so the main calculation can work with them. Lower them for systems with elements in traces, so the start does not disturb the element totals.
     - Any small positive amount.
   * - ``pa_DFYs``
     - 1e-6
     - Amount (moles) at which a pure phase is added back when the solver decides it was wrongly left out. It is never more than the bulk composition can supply. Rarely needs changing.
     - Any small positive amount.
   * - ``pa_PLLG``
     - 30000
     - Checks whether the calculation is drifting off in the wrong direction and stops it if so. 1 to 1000 is the useful range for this check. Lower it to catch failing calculations earlier; 0 if the check stops calculations that would succeed.
     - 0 = no check; 30000 or more also allows full diagnostic tracing.
   * - ``pa_PSM``
     - 1
     - How many diagnostic messages are written. Use 2 or 3 when investigating a failed calculation; 0 for large batch runs.
     - 0 = none (no log file); 2 = also warnings; 3 = detailed trace.
   * - ``pa_DNS``
     - 12.05
     - Standard density of surface sites (per nm²), used for the activity of surface species in sorption models. Sorption models that use a different standard site density.
     - Any positive number.
   * - ``pa_IEPS``
     - 0.001
     - How precisely the surface terms of sorption models are computed. Tighten it for sorption systems that converge poorly.
     - From 0.01 to 1e-6.
   * - ``pa_DKIN``
     - 1e-10
     - Tolerance on species amounts that are held between an upper and a lower limit (for example in kinetic or metastability calculations). Kinetic or metastability calculations with very small allowed ranges.
     - Any small positive amount.

Settings added in GEMS3K 5.0, which older project files do not contain:

.. list-table::
   :header-rows: 1
   :widths: 17 8 45 30

   * - Setting
     - Default
     - What it does and when it is useful
     - Other values
   * - ``pa_ColdRetryNudges``
     - 4
     - If a calculation from scratch (AIA) fails, retry up to this many times with a microscopically changed composition, then finish at the exact one. Leave on. Systems where AIA fails occasionally.
     - 0 = no retry.
   * - ``pa_IpmStallWindow``
     - 30
     - Stops once the energy and composition have clearly stopped changing for this many iterations. Leave on. Saves iterations on systems that converge slowly at the end.
     - 0 = off.
   * - ``pa_IpmAugmentedKKT``
     - 2
     - How each step's equations are solved. 2 is the most accurate. Leave at 2; other values only to compare with older results.
     - 0 = the original method; 1 = an intermediate method.
   * - ``pa_MbReproject``
     - 1 (on)
     - A final touch-up that makes the element totals add up exactly. Also used after phase selection, where a partly successful touch-up is now kept because it only seeds the next iterations. Leave on. Turn off only to see the raw result.
     - 0 = off.
   * - ``pa_PSTALL``
     - 1 (on)
     - Lets the mass-balance step give up early when it stops improving. Leave on. Saves time on systems where that step stalls.
     - 0 = off.
   * - ``pa_MbClassRule``
     - 0 (off)
     - Checks the totals of major elements in absolute terms and of trace elements in relative terms; the value is the trace/major ratio. Systems with trace elements at very low amounts.
     - A positive ratio = on.
   * - ``pa_FilloutBudget``
     - 0 (off)
     - Limits how much the starting guess may disturb the element totals, as a fraction of each element's amount. Systems with elements present only in traces.
     - A positive fraction = on.
   * - ``pa_DeterminacyWarn``
     - 0.01
     - Warns when the amount of a phase in the answer is fixed by the energy only to worse than this relative uncertainty, i.e. when small amounts are not reliable. Leave on; it tells you when reported small amounts should not be relied on.
     - 0 = no warning.
   * - ``pa_StabTPD``
     - 1 (on)
     - Reports, in the diagnostic trace only, whether a mixed phase left out of the answer should have been present. Does not change the answer. Investigating whether a solid solution or melt was wrongly left out.
     - 0 = off.

These names are still accepted in project files, so existing files load unchanged, but they no longer have
any effect: ``pa_LpDualFillout``, ``pa_OptimaLSWindow``, ``pa_IpmLoopTweaks``, ``pa_OptimaFDDiagFloor``,
``pa_MbPivotSplit``, ``pa_OptimaPhaseCompaction``, ``pa_OptimaReadmitSeed``, ``pa_MbTrendPhaseDecay``,
``pa_OptimaMaxStepRatio``, ``pa_GAR``, ``pa_GAH``.


Optima solver settings
~~~~~~~~~~~~~~~~~~~~~~

The ones most often adjusted:

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Setting
     - Default
     - What it does
   * - ``pa_OptimaTol``
     - 1e-8
     - How precisely Optima must converge. Raise it slightly to trade accuracy for speed on large systems.
   * - ``pa_OptimaMaxSeconds``
     - 0 (off)
     - Time limit in seconds for one calculation including retries. Useful in long batch or transport runs.
   * - ``pa_OptimaDimReduce``
     - 0 (auto)
     - First solves with only the species that matter, then adds the rest. Automatic means on from 200 species.
   * - ``pa_OptimaColdRetry``
     - 2
     - If a warm-started calculation (``sop``) fails, start again from scratch.
   * - ``pa_OptimaFinish``
     - 1 (on)
     - If Optima stops just short of converging, finish using the phases it found.
   * - ``pa_OptimaAcceptRepair``
     - 0 (off)
     - Corrects a small error in the element totals before accepting a nearly converged answer. Helps systems with very little water.
   * - ``pa_OptimaZeroAbsent``
     - 2
     - How absent species are reported. Use 1 for exact zeros for absent phases.

The full list, including the fine-tuning settings, is in the GEMS3K documentation.


Known limitations
-----------------

* **Optima is slower per calculation** than the GEMS3K IPM solver, often by one to two orders of magnitude on
  large systems. Use it where it improves the answer, not by default.
* **Some strongly non-ideal systems** can still be difficult to solve.
* **Very little water.** Close to the point where water runs out, ``aop`` and ``sop`` can fail. Use ``aia`` or
  ``hop`` there. Slightly further from that point, they work if the project sets ``pa_OptimaAcceptRepair = 1`` and
  ``pa_GAS = 0.002``.
* **Redox not fixed by the system.** Without a redox couple the energy and amounts are well defined, but
  Eh is not, and it can differ between runs and modes. Compare solvers on energy and phases, not Eh.
* **No free water solution** (very dry systems): pH, Eh and ionic strength are not meaningful and can
  differ widely between solvers. The phases are still correct.
* **Trace elements:** their mass balance can be relatively off in the Optima modes while the energy and
  the main elements are exact.
* **GEMS3K IPM solver iteration counts** vary on many systems with tiny input changes. The answer does not.
* **Trace phases:** the GEMS3K IPM solver can leave out a slightly supersaturated phase that Optima finds.
* **Warm-started IPM (``sia``)** along a sequence can keep phases that are no longer stable when a phase
  boundary is crossed.
* **Large systems** can cost tens of seconds per point in the Optima modes, and some fail where the GEMS3K IPM
  solver succeeds.
