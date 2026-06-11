"""Sphinx configuration for the tcren documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

project = "tcren"
author = "Antigenomics"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "nbsphinx",
]

autosummary_generate = False
napoleon_use_ivar = True  # render dataclass "Attributes:" inline (avoids duplicate objects)
autodoc_member_order = "bysource"
autodoc_typehints = "description"
nbsphinx_execute = "never"

# Heavy / optional dependencies mocked at doc-build time.
autodoc_mock_imports = ["arda", "scipy", "Bio", "matplotlib"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
