Solver Settings
===============

This page lists the solver settings that a user can change, and the known limits of the solvers.
Most projects need no change. Read the :doc:`solver_guide` first.


.. _project-file-settings:

Settings in project files
-------------------------

Both solvers read their settings from the ``-ipm`` file of the project (``pa_...`` entries).
They read the file one time, when xGEMS loads the project. xGEMS has no call that changes a setting later.
To change a setting, edit the file and make a new engine.

A setting that is not in the file keeps its default. Files from earlier versions do not contain the newer settings.
Add them when you need another value.

**JSON project.** Put the setting in the ``"ipm"`` block of ``...-ipm.json``:

.. code-block:: json

   "pa_PLLG": 10000,
   "pa_OptimaTol": 1e-09,
   "pa_OptimaMaxSeconds": 60,

**Text project.** Put the setting in ``...-ipm.dat`` after the ``<END_DIM>`` line, for example ``<pa_OptimaTol>  1e-09``.
xGEMS does not read a setting before ``<END_DIM>``.

Write each setting one time in the file. If a name occurs twice, xGEMS uses the later value.
Some settings, such as ``pa_PE``, ``pa_DHB`` and ``pa_DG``, apply to both solvers.

For suggested settings in common cases (transport, grids, dry systems and others), see :ref:`suggested-settings`.


GEMS3K IPM solver settings
--------------------------

Most project files contain these settings.

.. list-table::
   :header-rows: 1
   :widths: 17 8 45 30

   * - Setting
     - Default
     - What it does and when to change it
     - Other values
   * - ``pa_DK``
     - 1e-6
     - Sets the convergence precision. A smaller value gives a more precise answer and more iterations. Use 1e-7 for solid solutions with almost equal end-members.
     - Any small positive number.
   * - ``pa_DHB``
     - 1e-13
     - Sets how closely the amount of each element must agree with the input, as a fraction of the input. Loosen it if a calculation fails only on this check. Tighten it for trace elements.
     - Usually 1e-9 (looser) to 1e-15 (stricter).
   * - ``pa_DT``
     - 0
     - Sets how the check above is applied. 0 applies the same relative check to each element. Change it for systems with major and trace elements.
     - -6 or less = for major elements, use a fixed amount.
   * - ``pa_IIM``
     - 7000
     - Sets the maximum number of iterations of the main calculation. Raise it for large systems that run out of iterations.
     - Up to 32767. A larger value becomes 32767, with a warning.
   * - ``pa_DP``
     - 130
     - Sets the maximum number of iterations of the mass-balance step. Raise it if this step runs out of iterations.
     - Any positive number.
   * - ``pa_DW``
     - 1 (on)
     - Treats a run out of iterations as an error. Switch it off only to look at a result that GEMS3K rejects.
     - 0 = no error.
   * - ``pa_PE``
     - 1 (on)
     - Requires an electrically neutral answer. Keep it on for systems with ions. Switch it off only for systems with no charged species.
     - 0 = off.
   * - ``pa_PC``
     - 2
     - Selects the method that finds the present phases. 2 also adds phases that the solver left out by mistake. Keep 2.
     - 1 = the old method, for comparison with old results.
   * - ``pa_DF``
     - 0.01
     - Sets how clearly a missing phase must be stable before the solver adds it. Lower it if phases are missing. Raise it if phases appear and disappear along a sweep.
     - Any small positive number. A smaller value adds phases more easily.
   * - ``pa_DFM``
     - 0.01
     - Sets how clearly a present phase must be unstable before the solver removes it. Change it together with ``pa_DF``.
     - Any small positive number.
   * - ``pa_PD``
     - 2
     - Sets how often the solver calculates the activity coefficients. 2 = at each iteration. Keep 2. A lower value can save time for almost ideal systems.
     - 0 = one time at the start. 1 = only in the mass-balance step. 3 = only in the main calculation.
   * - ``pa_AG``
     - 1
     - Damps the non-ideal corrections between iterations. Lower it if a strongly non-ideal system oscillates.
     - -1 to 1.
   * - ``pa_DGC``
     - 0
     - A second damping setting that works with ``pa_AG``. Use 0.001 for sorption systems.
     - -1 to 1.
   * - ``pa_PRD``
     - -5
     - Removes small leftover amounts of species at the end. -5 removes amounts below 1e-5. Use a more negative value if real trace amounts disappear. 0 shows the raw result.
     - 0 = off. -6 or less = remove only smaller amounts.
   * - ``pa_GAS``
     - 0.001
     - Sets how strict the final clean-up is. A larger value keeps more trace amounts. Use 0.002 (with ``pa_OptimaAcceptRepair``) for systems with very little water.
     - Any small positive number.
   * - ``pa_DS``
     - 1e-20
     - Sets the smallest phase amount (mol) that GEMS3K reports as present. Raise it to hide trace phases.
     - Any small positive number.
   * - ``pa_XwMin``
     - 1e-13
     - Removes the aqueous solution if the amount of water (mol) is below this value. Use it for systems that dry out.
     - Any small positive number.
   * - ``pa_ScMin``
     - 1e-13
     - Removes the sorption phase if the amount (mol) of its solid is below this value. Use it if the sorbent dissolves almost completely.
     - Any small positive number.
   * - ``pa_PhMin``
     - 1e-20
     - Removes a non-aqueous solution phase with all its species if its amount (mol) is below this value. Use it for many solid solutions or melts in trace amounts.
     - Any small positive number.
   * - ``pa_DcMin``
     - 1e-33
     - Removes a species in a mixed phase if its amount (mol) is below this value. Seldom needs a change.
     - Any small positive number.
   * - ``pa_DB``
     - 1e-17
     - Sets the smallest element amount (mol) in the bulk composition (without charge). Use it for elements at extremely low amounts.
     - Any small positive number.
   * - ``pa_ICmin``
     - 1e-5
     - Sets the ionic strength (molal) below which the activity coefficients of aqueous species are 1. Change it for very dilute waters.
     - Any small positive number.
   * - ``pa_DG``
     - 1000
     - Scales the system internally to this total amount (mol). GEMS3K reports the original amounts. Switch the scaling off only to compare with unscaled results.
     - Below 1e-4 = no scaling.
   * - ``pa_EPS``
     - 1e-10
     - Sets the precision of the automatic start point (aia). Loosen it if the start point is not found. Tighten it for trace elements.
     - 1e-6 (looser) to 1e-14 (stricter).
   * - ``pa_DFYw``, ``pa_DFYaq``, ``pa_DFYid``, ``pa_DFYr``, ``pa_DFYh``, ``pa_DFYc``
     - 1e-5
     - Set the small amounts (mol) that the solver gives to water, aqueous species, species of ideal and non-ideal mixed phases, and pure phases that are zero in the start point. Lower them for trace elements.
     - Any small positive number.
   * - ``pa_DFYs``
     - 1e-6
     - Sets the amount (mol) of a pure phase that the solver adds back after it left the phase out by mistake. The bulk composition limits the amount. Seldom needs a change.
     - Any small positive number.
   * - ``pa_PLLG``
     - 30000
     - Stops a calculation that drifts in the wrong direction. 1 to 1000 is the useful range. Lower it to catch failures sooner. Use 0 if the check stops good calculations.
     - 0 = no check. 30000 or more also gives full diagnostic tracing.
   * - ``pa_PSM``
     - 1
     - Sets how many diagnostic messages GEMS3K writes. Use 2 or 3 to find the cause of a failure. Use 0 for large batch runs.
     - 0 = none (no log file). 2 = also warnings. 3 = detailed trace.
   * - ``pa_DNS``
     - 12.05
     - Sets the standard density of surface sites (per nm²) for the activity of surface species. Change it if a sorption model uses another density.
     - Any positive number.
   * - ``pa_IEPS``
     - 0.001
     - Sets the precision of the surface terms of sorption models. Tighten it if a sorption system converges badly.
     - 0.01 to 1e-6.
   * - ``pa_DKIN``
     - 1e-10
     - Sets the tolerance for species amounts that have an upper and a lower limit (kinetic or metastability calculations). Use it if the allowed range is very small.
     - Any small positive number.

GEMS3K 5 added these settings. Older project files do not contain them.

.. list-table::
   :header-rows: 1
   :widths: 17 8 45 30

   * - Setting
     - Default
     - What it does and when to change it
     - Other values
   * - ``pa_DeterminacyWarn``
     - 0.01
     - Writes a warning if the energy fixes the amount of a phase with a relative uncertainty above this value. Such small amounts are not reliable. Keep it on.
     - 0 = no warning.
   * - ``pa_FailRescue``
     - 2 (on)
     - Both solvers. Tries a failed or BAD call again in steps. Keep 2 for transport codes. A call that no step can save stops after ``pa_RescueBudget`` (100000) iterations. See :ref:`batch-and-transport`.
     - 0 = off. 1 = steps 1 to 5 only.


Optima solver settings
----------------------

A cold calculation (``aop``) first tries a quick solve. If the quick solve fails, the full solve runs.
Optima stops a solve that makes no progress and tries again. A failed warm calculation (``sop``) starts again cold.

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Setting
     - Default
     - What it does and when to change it
   * - ``pa_OptimaTol``
     - 1e-8
     - Sets the convergence precision. A larger value gives less accuracy and more speed on large systems.
   * - ``pa_OptimaMaxSeconds``
     - 0 (off)
     - Sets the time limit (seconds) for one Optima solve. The results then depend on the speed of the computer. The limit starts again for each try and each rescue step.
   * - ``pa_OptimaCallCap``
     - 0 (no cap)
     - Sets the largest number of iterations for one ``aop`` or ``sop`` calculation. After the cap, GEMS3K repeats the calculation with amounts that it changes by a very small value. -1 = the cap depends on the size of the system.
   * - ``pa_OptimaDimReduce``
     - 0 (auto)
     - Solves first with only the important species. Then it adds the other species. Auto = on from 200 species. A negative number = off.
   * - ``pa_OptimaColdRetry``
     - 2
     - Starts again cold after a failed warm calculation (``sop``, ``shp``). 2 = stop the first try early. 1 = start again only after the full budget. 0 = no new try.
   * - ``pa_OptimaAcceptRepair``
     - 0 (off)
     - Corrects a small error in the element totals before it accepts a nearly converged result. It helps systems with very little water.
   * - ``pa_OptimaZeroAbsent``
     - 2
     - Sets how the result shows absent species. 1 = exact zeros for absent phases.
   * - ``pa_OptimaRedoxFree``
     - 1
     - Writes a log message if the bulk composition does not fix the redox level. 2 = also aim at pe 4 when you give no pH or Eh target. 0 = off.
   * - ``pa_OptimaDryOut``
     - 0 (off)
     - Solves the solids without the water phase when the water is gone. Use it only if the system uses all the water by design (for example cement hydration).

The GEMS3K documentation gives more detail:

* The `options guide <https://github.com/gemshub/GEMS3K/blob/develop_optima2/Docs/SOLVER-OPTIONS-GUIDE.md>`_ lists all options.
* The `solver path guide <https://github.com/gemshub/GEMS3K/blob/develop_optima2/Docs/SOLVER-PATH-GRAPH.md>`_ shows the order of each step and try.
* The `glossary <https://github.com/gemshub/GEMS3K/blob/develop_optima2/Docs/GLOSSARY.md>`_ defines the technical words.

Settings for tests are on the page :doc:`expert_settings`.


Known limits
------------

* **Optima is slower for each calculation.** It is slower than the GEMS3K IPM solver, and much slower on large systems with many solid solutions. Use Optima where it gives a better result. Do not use it as the default.
* **Strongly non-ideal systems** can be difficult to solve, also with Optima.
* **Very little water.** Near the point where the water is gone, the modes can give different trace phases. Optima can need many more iterations.
  If the system uses all the water by design, set ``pa_OptimaDryOut = 1``. The result has no water phase, and the pH is 0.
  If ``aop`` or ``sop`` fails, set ``pa_OptimaAcceptRepair = 1`` and ``pa_GAS = 0.002``, or use ``aia`` or ``hop``.
* **Redox level not fixed.** A system without a redox couple has defined energy and amounts, but no defined Eh. Eh can differ between runs and modes.
  Compare the energy and the phases, not the Eh. For a defined Eh, set ``pa_OptimaRedoxFree = 2`` (it aims at pe 4) or give an Eh target.
* **Small input changes** can change the pH and the Eh. The phases and the energy do not change.
* **No free water solution** (very dry systems). The pH, the Eh and the ionic strength have no meaning. They can differ much between solvers. The phases are correct.
* **Trace elements.** In a few cases the relative mass balance of trace elements is not good in the Optima modes. The energy and the main elements are exact.
* **Trace elements below the species floors.** The floor amounts of the species of a trace element can add up to more than the amount of the element.
  An Optima call can then stay BAD, and the rescue cannot solve it. Use ``aia`` or ``hop``, or raise ``min_amount`` (see :doc:`advanced_usage`).
* **Solid solutions near a miscibility gap, and twin phases,** need many more iterations. The labels of twin phases can change between modes. The energy is the same.
* **GEMS3K IPM solver iteration counts** change on many systems when the input changes by a very small amount. The result does not change.
* **Trace phases.** The GEMS3K IPM solver can leave out a phase that is a little supersaturated. If GEMS3K has Optima, it finishes such results with Optima (``pa_NativeOptimaFinish``).
* **A warm ``sia`` sequence** can keep phases that are no longer stable after a phase boundary. Use ``aop`` or ``sop``.
* **Large systems** can need a long time for each point in the Optima modes.
* **Rescued calls cost more.** A call that the rescue saves can use several times the iterations of a normal call.
  A call that no step can save uses the most. One such call can take more time than all other points of a grid.
  ``pa_RescueBudget`` limits the iterations of such a call. ``pa_OptimaMaxSeconds`` does not limit the total time,
  because it starts again for each step. For grids where an error is acceptable, see :ref:`suggested-settings`.
