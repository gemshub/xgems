Expert Settings
===============

.. warning::
   These settings are for tests of GEMS3K. Do not change them. Use the default value.

Each setting switches a correctness fix or a speed tuning on or off, or sets an internal limit.
The settings are separate so that a test can check each mechanism alone.
A project file keeps a setting only if its value is not the default.
The same list is in the GEMS3K
`advanced settings <https://github.com/gemshub/GEMS3K/blob/develop_optima2/Docs/SOLVER-ADVANCED-SETTINGS.md>`_.

The settings that a user can change are on the page :doc:`solver_settings`.
See :ref:`project-file-settings` for how to put a setting in a project file.


Finish and acceptance of a result (Optima)
------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_OptimaFinish``
     - 1 (on)
     - Optima can stop just short of convergence. Then a final step finishes the calculation on the phases that Optima found.
   * - ``pa_OptimaFinishFirst``
     - 1 (on)
     - After a stuck calculation fails, the solver runs the finishing step on the stop point first. It keeps that result if no phase tends to form.
   * - ``pa_OptimaFinishStall``
     - 15 (on)
     - The finishing step stops after this number of steps without progress.
   * - ``pa_OptimaLevelFinish``
     - 1 (on)
     - The leveled start of a cold calculation can hold mixed phases. Then the finishing step runs on that start first.
   * - ``pa_OptimaTpdAccept``
     - 1e-6 (on)
     - Accepts a result that is not fully converged. A direct check must show that no missing phase must be present (within this tolerance).
   * - ``pa_OptimaTpdVeto``
     - 1 (on)
     - The direct check of ``pa_OptimaTpdAccept`` cannot overrule the phase stability check. This is for a mixed phase that the direct check cannot test, for example a water phase that is gone.
   * - ``pa_OptimaCollapseCheck``
     - 1 (on)
     - The check of a result judges a phase that disappeared by its pure species. A phase that must come back is then not accepted as absent. A pure solid below the presence limit (``pa_DS``) is judged by its chemical potential. A water phase that has almost vanished is judged as absent.
   * - ``pa_OptimaTwinReadmit``
     - 1 (on)
     - The solver can switch off a mixed phase because a copy of it took over. After the calculation, the solver tries the phase again if the phase stability check says that it must form.
   * - ``pa_OptimaTwinLabel``
     - 1 (on)
     - Two copies of one solid solution can be present at two different compositions. After an Optima calculation, each composition goes into the copy whose major end-members (class M) it holds most. This is done only when the two copies have the same species, standard energies, limits and mixing parameters. G does not change. 0 = off.
   * - ``pa_OptimaTraceRespread``
     - 1 (on)
     - After the calculation, an element that only the floor amounts of its species hold gets its real amount back.
   * - ``pa_OptimaWaterKeep``
     - 2 (on)
     - The value 1 keeps a converged result when the new try with more water fails.


Line search and curvature (Optima)
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_OptimaLineSearch``
     - 1.5
     - A step can make the result worse by more than this factor. Then the solver tries a shorter step.
   * - ``pa_OptimaLSStallEscape``
     - 10 (on)
     - Shorter steps can fail this number of times in a row. Then the solver takes one full step.
   * - ``pa_OptimaLSRejectWorse``
     - 0 (off)
     - If a shorter step does not help, the solver uses the full step.
   * - ``pa_OptimaLSWarmTrials``
     - 5 (on)
     - In a warm calculation (SOP, SHP), the line search tries at most this number of step lengths.
   * - ``pa_PhaseHessianFloor``
     - 0.01 (on)
     - Keeps the curvature of mixed phases in a good range.
   * - ``pa_OptimaFDHessian``
     - 1 (on)
     - Calculates a more exact curvature for non-ideal phases. This costs more time.
   * - ``pa_OptimaFDHessianDelay``
     - 0
     - Uses the cheaper approximate curvature for this number of first iterations.
   * - ``pa_OptimaMoleFracHessian``
     - 0 (off)
     - Uses a more exact curvature for ideal mixed phases.
   * - ``pa_OptimaFloorCurvature``
     - 3 (on)
     - Uses the exact curvature of the logarithm term for each species in a mixed phase. The numerical estimate is too low for very small amounts.
   * - ``pa_OptimaFDPureSkip``
     - 1 (on)
     - The solver sets the curvature of a pure phase to exactly zero. It does not calculate an activity.
   * - ``pa_OptimaPartialPivLU``
     - 100 (on)
     - Optima uses a faster method for its linear equations when there are at least this number of equations.


First solve and cold start (Optima)
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_OptimaFastFirst``
     - 60 (on)
     - A cold calculation (AOP) first tries a quick solve with full steps. If the quick solve does not end well, the normal calculation runs from the same start.
   * - ``pa_OptimaFastFrozen``
     - 3 (on)
     - The quick first try stops when this number of steps in a row change almost nothing.
   * - ``pa_OptimaFastDiverge``
     - 10 (on)
     - The quick first try stops when its energy jumps far above the start and its error does not improve by this factor.
   * - ``pa_OptimaFastRises``
     - 10 (on)
     - The quick first try stops after this number of steps where the energy rises by more than 1 %. Then the full calculation runs.
   * - ``pa_OptimaFastCollapse``
     - 1000 (on)
     - The quick first try stops at once when its energy is more than this number of times its start value.
   * - ``pa_OptimaFastAcceptCut``
     - 1 (on)
     - The quick first try can use all its iterations and still end OK.
   * - ``pa_OptimaDimReduceTol``
     - 10
     - Sets how many species the first guess of "important species" includes. A larger value includes more species.
   * - ``pa_OptimaPreSolveFirstIters``
     - 6000
     - Sets the iteration limit for each pass of the reduced first solve. The solver stops a pass that makes no progress.
   * - ``pa_OptimaPreSolveSeed``
     - 1 (on)
     - The first pass of the reduced first solve starts from the simple energy-minimizing estimate. If a pass needs more than 300 steps, the solver does not use the reduced first solve for that calculation.
   * - ``pa_OptimaEmptyPreSolveSkip``
     - 1 (on)
     - The reduced first solve can find no usable set of species. Then the new tries of that calculation do not use it.
   * - ``pa_OptimaColdLevel``
     - 1 (on)
     - If the system does not use the reduced first solve, a cold calculation starts from a simple energy-minimizing estimate. Two copies of one solid solution start at different compositions.
   * - ``pa_OptimaLevelStart``
     - 1 (on)
     - For systems with only solids and other mixed condensed phases, a cold calculation first finds which mixed phases form and their compositions. It then starts from that result.
   * - ``pa_OptimaCgSeed``
     - 1e-6 (on)
     - Gives a better start to a new cold try. It uses the phases that are likely to be stable.
   * - ``pa_OptimaEarlyStabilityAt``
     - 0 (auto)
     - Lets Optima find early that a phase must disappear. Without it, the phase becomes smaller slowly.
   * - ``pa_OptimaDcFloor``
     - 0 (auto)
     - Sets the smallest amount that Optima can represent for an absent species.


Stop rules (Optima)
-------------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_OptimaStallWindow``
     - 200 (on)
     - Stops a calculation that does not improve for this number of iterations. A new try can then start sooner.
   * - ``pa_OptimaStallNetFactor``
     - 5 (on)
     - After a stall, Optima solves one more time with the stall check off. This last solve gets at most this number times the iterations of the solve that stalled, and at least 2000. 0 = the normal limit.
   * - ``pa_IpmTwinSplit``
     - 1e-6 (on)
     - After a warm GEMS3K IPM calculation, two copies of one solid solution can keep one composition where they should split into two phases. The solver then tests that composition. If a split is lower by more than this value (RT per mole), the calculation is solved again from a cold start. 0 = off.
   * - ``pa_OptimaFrozenCut``
     - 20 (on)
     - Stops a calculation when its energy and error do not change for this number of iterations. The finishing step and the new tries then take over.
   * - ``pa_OptimaTailStop``
     - 100 (on)
     - Stops a calculation when its energy does not change for this number of iterations and the error is only on trace amounts. The solver checks the phase stability at the stop point. If the check passes, the solver takes the point as the result.
   * - ``pa_OptimaTailStopWarm``
     - 10 (on)
     - The same as ``pa_OptimaTailStop``, for a warm calculation (SOP, SHP).
   * - ``pa_OptimaTailDryOnly``
     - 2 (on)
     - Sets where the tail stop acts.
   * - ``pa_OptimaSnapLeaving``
     - 2 (on)
     - In a warm restart, a species that goes to zero can limit the step many times in a row. After this number of steps, the solver puts the species on its floor.
   * - ``pa_OptimaCapNudges``
     - 3 (only with ``pa_OptimaCallCap``)
     - Sets how many times the solver repeats a calculation that reached its cap. Each repeat changes the amounts by a very small value.


Water (Optima)
--------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_OptimaEmptySolvent``
     - 1 (on)
     - If the water in a cold calculation is at or below its minimum amount, the solver adds water. This is the case when the dissolved amounts are larger than the water.
   * - ``pa_OptimaDryRetry``
     - 0 (off)
     - A calculation can end with almost no liquid water, while the bulk can form much more. Then the solver solves again from a start that first puts hydrogen and oxygen into water.


Both solvers
------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_TpdStrongBasins``
     - 4 (on)
     - The thorough check of whether a mixed phase can form tries each start composition one time.
   * - ``pa_LogBarrierTau``
     - 1e-16
     - A very small numerical push that keeps pure minerals away from exactly zero.
   * - ``pa_StabTPD``
     - 1 (on)
     - Reports, in the diagnostic trace only, that a mixed phase must be in the result when it is not.
   * - ``pa_FilloutBudget``
     - 0 (off)
     - Limits how much the start guess can change the element totals. The limit is a fraction of the amount of each element.
   * - ``pa_RescueWalkSteps``
     - 8
     - Sets the number of equal steps in the walk from the previous result to the new composition.
   * - ``pa_RescueBudget``
     - 100000 (on)
     - The rescue (``pa_FailRescue``) skips steps 1 to 4 after the call has used this number of iterations. The count includes the first attempt and both solvers. Step 5, the other solver, still runs one time. 0 = no limit.


GEMS3K IPM solver
-----------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Setting
     - Default
     - What it does
   * - ``pa_ColdRetryNudges``
     - 4 (on)
     - If a cold calculation (AIA) fails, the solver tries again up to this number of times. Each try uses a composition that is changed by a very small value. The solver then finishes at the exact composition.
   * - ``pa_IpmStallWindow``
     - 30 (on)
     - Stops when the energy and the composition do not change for this number of iterations.
   * - ``pa_IpmNoProgress``
     - 500 (on)
     - Needs a build with Optima. The solver stops if the energy does not improve for this number of iterations.
   * - ``pa_IpmAugmentedKKT``
     - 2
     - Sets how the solver solves the equations of each step.
   * - ``pa_IpmLevelStart``
     - 2 (on)
     - A cold calculation (AIA) on a system without aqueous, gas or sorption phases starts from the leveled solution. This is a linear calculation with the mixed phases at their best compositions.
   * - ``pa_IpmQrSortPivot``
     - 100 (on)
     - A step with at least this number of active species uses a faster and more accurate method for its equations.
   * - ``pa_MbReproject``
     - 1 (on)
     - A final correction that makes the element totals agree exactly.
   * - ``pa_PSTALL``
     - 1 (on)
     - Lets the mass-balance step stop early when it does not improve.
   * - ``pa_MbClassRule``
     - 0 (off)
     - Checks the totals of major elements in absolute terms and the totals of trace elements in relative terms. The value is the ratio of trace to major.
   * - ``pa_NativeOptimaFinish``
     - 500 (on)
     - Needs a build with Optima. A solid can be left out, or be in a trace amount, and still tend to form. Then Optima finishes the result. The solver keeps the new result only if it is better. Optima also finishes a calculation that reaches its iteration limit or fails to solve its equations, from the point where it stopped.


Names with no effect
--------------------

xGEMS accepts these names in project files, so old files load unchanged. They have no effect:
``pa_LpDualFillout``, ``pa_OptimaLSWindow``, ``pa_IpmLoopTweaks``, ``pa_OptimaFDDiagFloor``,
``pa_MbPivotSplit``, ``pa_OptimaPhaseCompaction``, ``pa_OptimaReadmitSeed``, ``pa_MbTrendPhaseDecay``,
``pa_OptimaMaxStepRatio``, ``pa_GAR``, ``pa_GAH``.
