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
author = "ISALGO laboratory"

# Single-source the version so the published docs never drift (the old hardcoded "0.1.0" was left
# stale while the package moved to 2.1.x). Read it without importing tcren — its heavy deps are
# only mocked during the autodoc pass, not at conf-import time — so parse pyproject.toml, which is
# where the version now lives; `tcren.__version__` reads it back out of the installed metadata, and
# a checkout that was never installed reports 0.0.0+unknown, which is not what the docs should say.
import re as _re  # noqa: E402

_pyproject = (Path(__file__).resolve().parent.parent / "pyproject.toml").read_text()
release = _re.search(r'^version\s*=\s*["\']([^"\']+)["\']', _pyproject, _re.M).group(1)
version = ".".join(release.split(".")[:2])

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
autodoc_mock_imports = [
    "arda", "scipy", "Bio", "matplotlib",
    "py3Dmol", "openmm", "pdbfixer", "promod3", "ost",  # optional viz / refinement engines
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
