# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import site

# xgems.PyxGEMS is a compiled pybind11 extension shipped with a PEP 561 .pyi
# stub (for IDE/type-checker support). Sphinx >= 9's autodoc prefers that
# stub over the real module for native extensions, but the stub has no
# docstrings (just `...` bodies), so every member ends up undocumented.
# This tells autodoc to import the real .so and use its actual docstrings.
os.environ['SPHINX_AUTODOC_IGNORE_NATIVE_MODULE_TYPE_STUBS'] = '1'

sys.path.insert(0, os.path.abspath('../../..'))
sys.path.insert(0, os.path.abspath('../../../python'))
sys.path.insert(0, site.getsitepackages()[0])

project = 'xGEMS'
copyright = '2025, GEMS Team'
author = 'Allan Leal, Dmitrii Kulik, G.D. Miron'
release = '1.0.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode"]

autodoc_default_options = {
    'exclude-members': 'pybind11_detail_function_record_v1_system_libstdcpp_gxx_abi_1xxx_use_cxx11_abi_1',
}
breathe_projects = {
    "xGEMS": "../../build/doxygen/xml/"
}
breathe_default_project = "xGEMS"



# Enable automatic Python documentation extraction
# autodoc_mock_imports = ["xgems.PyxGEMS"]

#templates_path = ['_templates']
#exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
#html_static_path = ['_static']

html_theme_options = {
    "collapse_navigation": False,
    "sticky_navigation": True,
}

import subprocess
subprocess.call('cd ../../.. ; doxygen', shell=True)

html_extra_path = ['../../html']
