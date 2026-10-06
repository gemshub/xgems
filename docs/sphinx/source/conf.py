# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import re
import subprocess
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

def _is_prerelease():
    """True when the docs are built from a branch other than the main one, so that the
    prerelease install instructions (``only:: prerelease``) appear only there."""
    forced = os.environ.get('XGEMS_DOCS_PRERELEASE')
    if forced is not None:
        return forced == '1'
    name = os.environ.get('READTHEDOCS_VERSION_NAME')
    if name is not None:  # Read the Docs: "latest" is the default (main) branch
        return os.environ.get('READTHEDOCS_VERSION_TYPE') == 'branch' and name not in ('latest', 'stable')
    try:
        branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                                         stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return False
    return branch not in ('master', 'main', 'HEAD')


if _is_prerelease():
    tags.add('prerelease')  # noqa: F821  (``tags`` is provided by Sphinx)

project = 'xGEMS'
copyright = '2025, GEMS Team'
author = 'Allan Leal, Dmitrii Kulik, G.D. Miron'
release = '2.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode"]

autodoc_default_options = {
    'exclude-members': 'pybind11_detail_function_record_v1_system_libstdcpp_gxx_abi_1xxx_use_cxx11_abi_1',
}


def _drop_overload_signature(app, what, name, obj, options, lines):
    """pybind11 puts `name(*args, **kwargs)` first in the docstring of an overloaded function;
    reST reads the stars as unterminated emphasis, so remove that line."""
    if lines and re.match(r'^\w+\(\*args, \*\*kwargs\)$', lines[0].strip()):
        del lines[0]


def setup(app):
    app.connect('autodoc-process-docstring', _drop_overload_signature)


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

subprocess.call('cd ../../.. ; doxygen', shell=True)

# docs/html is the pre-built site that older branches publish; Read the Docs builds from source
if not os.environ.get('READTHEDOCS'):
    html_extra_path = ['../../html']
