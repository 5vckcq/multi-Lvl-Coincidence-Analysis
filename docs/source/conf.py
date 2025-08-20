# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../../src/mlca'))

import mlca


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'multi-Lvl-Coincidence Analysis'
copyright = '2025, Johannes Mierau'
author = 'Johannes Mierau'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.autosummary',
              'sphinx.ext.intersphinx', # Cross-references to other projects
              "sphinx.ext.coverage", # Collects doc coverage stats
              #"sphinx.ext.viewcode", # Links to highlighted source code (i.e. "[source]" button)
              'autoapi.extension',
              #"sphinx_automodapi.automodapi", # Automatically generates module documentation
              #"sphinx_automodapi.smart_resolver", # Helps resolving some imports
              "numpydoc", # Support for the NumPy docstring format
              ]
autosummary_generate = True
autoapi_dirs = ['../../src/mlca']
autodoc_mock_imports = ["pandas"]

templates_path = ['_templates']
exclude_patterns = ["config", "csv-files", "samples", "output"]

html_domain_indices = True

html_theme_options = {
    'nosidebar': False,
    'logo': "mlca-logo.png",
    "description" : "Causal-mechanistic modeling using Coincidence Analysis",
    'github_user': '5vckcq',
    'github_repo': 'multi-Lvl-Coincidence-Analysis',
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
