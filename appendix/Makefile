# Compile the tcren methods appendix. Requires a TeX distribution (latexmk + lualatex + bibtex);
# lualatex is needed for the OldStandard OpenType text+math fonts and is pinned in .latexmkrc.
# The workflow figure is a graphviz diagram (workflow.dot -> workflow.pdf); both are committed.
tcren.pdf: tcren.tex refs.bib workflow.pdf
	latexmk -lualatex -interaction=nonstopmode -halt-on-error tcren.tex

workflow.pdf: workflow.dot
	dot -Tpdf workflow.dot -o workflow.pdf

.PHONY: clean
clean:
	latexmk -C tcren.tex
