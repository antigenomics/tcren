# Minimal sphinx Makefile.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Regenerate the executed .ipynb of every marimo notebook, so nbsphinx (which runs with
# nbsphinx_execute = "never") has figures to render. Needs the structure data the notebook loads;
# the committed .ipynb is what a docs build without that data falls back on.
MARIMO_NOTEBOOKS = surface_topology

.PHONY: notebooks
notebooks:
	@for nb in $(MARIMO_NOTEBOOKS); do \
		echo "exporting $$nb"; \
		MPLBACKEND=Agg python -m marimo export ipynb --include-outputs --sort top-down \
			../notebooks/$$nb.py -o ../notebooks/$$nb.ipynb -f; \
	done

.PHONY: help Makefile notebooks

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
