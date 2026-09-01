# The descriptor whitepaper. `make` regenerates the catalogue layer, then builds the PDF.
TEX := descriptors.tex

all: descriptors.pdf

generated/catalogue.tex generated/reduction.pdf: ../src/tcren/descriptors/catalogue.py
	cd .. && python scripts/gen_appendix.py

descriptors.pdf: $(TEX) generated/catalogue.tex generated/reduction.pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(TEX)

check:
	cd .. && python scripts/gen_appendix.py --check

clean:
	latexmk -C $(TEX)

.PHONY: all check clean
