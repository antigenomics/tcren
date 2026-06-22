# Compile the tcren methods appendix. Requires a TeX distribution (latexmk + lualatex + bibtex);
# lualatex is needed for the OldStandard OpenType text+math fonts and is pinned in .latexmkrc.
tcren.pdf: tcren.tex refs.bib
	latexmk -lualatex -interaction=nonstopmode -halt-on-error tcren.tex

.PHONY: clean
clean:
	latexmk -C tcren.tex
