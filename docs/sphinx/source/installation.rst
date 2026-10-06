Installation
============

xGEMS is installed with `conda <https://docs.conda.io/en/latest/miniconda.html>`_ (Miniconda or Miniforge).
Create an environment and install the package from conda-forge:

.. code-block:: bash

   conda create -n xgems -c conda-forge xgems
   conda activate xgems

Check that it works:

.. code-block:: bash

   python -c "from xgems import ChemicalEngine; print('xGEMS is ready')"

To add xGEMS to an environment you already have, run ``conda install -c conda-forge xgems`` in it.
To build it from source, follow the instructions in the
`README <https://github.com/gemshub/xgems#readme>`_.

Examples of using xGEMS, as Jupyter notebooks, are in the
`xgems-jupyter <https://github.com/gemshub/xgems-jupyter>`_ repository.

.. only:: prerelease

   Install the prerelease
   ----------------------

   This documentation is for a **test version** of xGEMS that has the new Optima solver and is not
   released yet. You can try it without touching your normal xGEMS: it goes into its own environment.

   1. Create the environment. The command is long, but you only copy it:

      .. code-block:: bash

         conda create -n xgems-optima --override-channels --strict-channel-priority \
             -c https://gemshub.github.io/xgems/prerelease -c conda-forge xgems

   2. Switch to it:

      .. code-block:: bash

         conda activate xgems-optima

   3. Check that the new solver is there. This should print ``True``:

      .. code-block:: bash

         python -c "from xgems import ChemicalEngine; print(ChemicalEngine.builtWithOptima())"

   Then see the :doc:`solver_guide` to choose a solver.

   Good to know:

   * Always use a separate environment for the test version, as above, and not your normal one.
   * To remove it again: ``conda env remove -n xgems-optima``.
   * It is built for Python 3.12 to 3.14 on Linux, macOS (Apple silicon) and Windows. The macOS and Windows
     builds are still being checked.
   * If you were given a folder with the packages instead of the web address, use
     ``-c file:///path/to/folder`` in place of the web address (on Windows ``file:///C:/path/to/folder``).
   * The web address is provisional until the first test version is published. This page disappears from
     the released documentation.
