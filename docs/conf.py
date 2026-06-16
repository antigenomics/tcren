"""Sphinx configuration for the tcren documentation."""

import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("../src"))

# Tutorial notebooks live in the repo-root ``notebooks/`` directory; copy them into
# ``docs/notebooks/`` (gitignored) at build time so nbsphinx resolves their extracted
# images correctly (a symlink mis-resolves those paths).
_HERE = Path(__file__).resolve().parent
_NB_SRC = _HERE.parent / "notebooks"
_NB_DST = _HERE / "notebooks"
if _NB_SRC.is_dir():
    _NB_DST.mkdir(exist_ok=True)
    for _nb in _NB_SRC.glob("*.ipynb"):
        shutil.copy2(_nb, _NB_DST / _nb.name)

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
